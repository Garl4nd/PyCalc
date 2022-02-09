from abc import ABC,abstractmethod,abstractproperty
import itertools as it
import functools as ft
import cmath
import numpy as np
import scipy.special
from scipy.special import gamma as cgamma,erf,erfc,binom,beta,ellipkinc,ellipeinc
import unicodedata
from operator import mul,add,pow,truediv
module_prefix="" #module_prefix=module_prefix_optinos["direct"]  to make eval the inverse of repr

# TODO: komplexní koeficienty u monomů (místo genobj...) -> obecněji to mohou být obecné koeficienty (z GeneralObject)
# TODO: Řešení rovnic v závislosti na parametru - bylo by nějak těžké přesunout členy monomu ze slovníku do koeficientu? 
#       Nemělo by být, problém je v tom, že násobení a sčítání je definované přes monomy, takže by se to (asi) muselo ošetřit.
# TODO: Doladit Draw_Parse_Tree - opravit indexy a kolize unárních a binárních operátorů (je třeba tam udělat víc místa, např 2-#4)
# TODO: soustavy nelineárních rovnic
# TODO: řešení obecných rovnic aproximativními metodami
# TODO: krácení polynomů ve více proměnných - celkem challenge
# TODO: Ten __getitem__ se mi moc nelíbí, např. teď máme víc způsobů jak získat real a imag z komplexního čísla. Plus tam teď máme derivování. Aby to bylo čistší, tak by se mohlo přesunout z getitem do funkce v regparseru                       
# TODO: Počítání s nebodovýma funkcema: skládání, derivace, evaluace...jinak asi nepůjdou *parciální* diferenciální rovnice udělat pořádně
#Tím se ale dostáváme k násobení Jacobiánů...
# TODO: Přidat class derivace: class Deriv(GeneralObject), evaluables: GeneralObject , def eval(self,*dervars): return self.attached.deriv(*dervars), jednak kvůli dif. rovnicím, jednak třeba kvůli reprezentaci derivace dot productu
# TODO: Přidat funkce pro: transpozice, divergence, curl,matmul, cross
# Diskutabilní řádka: if all(isinstance(comp,ImplicitVector) for comp in values):
import math
def ignore(*ignorelist_str):
    def decorator(operator):
        
        @ft.wraps(operator)
        def func(x,y):  
            #print(x,y,type(y),type(y) in ignorelist,ignorelist[0])
            ignorelist=tuple(eval(el) for el in ignorelist_str)
            if any(isinstance(y,igntype) for igntype in ignorelist):
             #   #print("ignoring")
                return NotImplemented
            else:
                return operator(x,y)
        return func
    return decorator

def numpyconst(c):
    def const(*args,**kwargs):
        var=next(iter(kwargs))
        #print("npconst",np.shape(kwargs[var]),type(c))
        return c*np.ones(np.shape(kwargs[var]))
    return const

class GeneralObject(ABC):
    @abstractproperty
    def name(self):
        pass
    @abstractproperty
    def evaluables(self):
        return int,float
    @abstractproperty
    def is_evaluable(self):
        return False
    @property
    def tol(self):
        return 1.e-10

    def __init__(self,*attached_objs: "GeneralObject"):
        self.attached=attached_objs
        #print(self.attached)
        #self.tol=1.e-10
        super().__init__()
    
    @ignore("Equation")
    def __add__(self,other):
        
        
        if isinstance(other,(Monom,Polynomial,Quotient)):
            return  other+Monom(self)
        elif isinstance(other,Complex):
                return Complex(self+other.real,other.imag)
        elif isinstance(other,Vector):
            return Vector(*(self+comp for comp in other.values))
        else:
            #print([(str(el),repr(el)) for el in (Monom(self),Monom(other))])
            return Monom(self)+Monom(other)
        if self==other:
         #   #print("attached:",self.attached)
            new_obj=self.__class__(*self.attached)
            
          #  #print("returning ",new_obj)
            return new_obj
        return ObjSum(self,other)
    __radd__=__add__
    @ignore("Equation")
    def __sub__(self,other): 
        #print("subing",self,other)
        return self+(-1)*other
    def __rsub__(self,other):
        return -(self-other)
    
    @ignore("Equation")
    def __mul__(self,other):
        #print("str",self,other,type(self),type(other),isinstance(self,GeneralObject),isinstance(other,GeneralObject))
        #import genpolynom
        #if isinstance(other,(int,float)):
        #    return other*Monom(self)
        if isinstance(other,(int,float,Monom,Polynomial,Quotient)):
            #print("here",other,type(other))
            return other*Monom(self)
        if isinstance(other,Complex):
            return Complex(self*other.real,self*other.imag)
        if isinstance(other,Vector):
            return Vector(*(self*comp for comp in other))
        if isinstance(other,GeneralObject):
            return Monom(self)*Monom(other)
            
        return ObjMul(self,other)
    
    __rmul__=__mul__
    def __pow__(self,exponent):
        if isinstance(exponent,Vector):
            return exponent.__class__(*(self**comp for comp in exponent))
        if isinstance(exponent,(int,float)):
            return Monom(self)**exponent
        else:
            return ObjPow(self,exponent)
    def __rpow__(self,base):
#        if isinstance(base,(int,float))
        return ObjPow(base,self)
    @ignore("Equation")
    def __truediv__(self,other):
        if isinstance(other,Complex):
            return self*(1/other)
        if isinstance(other,Vector):
            return Vector(*(self/comp  for comp in other.values))

        return Polynomial(self)/other
    def __rtruediv__(self,other):
        if self==0:
            raise ValueError("Can't divide by zero!")
        else:
            
            return other/Monom(self)
    def __neg__(self):
        return (-1)*self
    @abstractmethod 
    def eval(self,*args,**kwargs):
        pass
    def get_vars(self):
        vars=set()
        for arg in self.attached:
            #print(self,arg,type(arg))
            if isinstance(arg,GeneralObject):
                
                vars.update(arg.get_vars())
        return vars
    def get_const_vars(self):
        vars=set()
        for arg in self.attached:
            #print(self,arg,type(arg))
            if isinstance(arg,GeneralObject):
                vars.update(arg.get_const_vars())
        return vars
        
    def contains_complex(self):
        for arg in self.attached:
            #print(self,arg,type(arg))
            if isinstance(arg,GeneralObject):
                if arg.contains_complex():
                    return True
        return False
    def contains_const(self):
        for arg in self.attached:
            #print(self,arg,type(arg))
            if isinstance(arg,GeneralObject):
                if arg.contains_const():
                    return True
        return False
    def deriv(self,dervar):
        pass
    def grad(self):
        vars=list(sorted(self.get_vars()))
        if vars:
            return Vector(*(self.deriv(var) for var in vars))
        else:
            return Monom(0)
    def laplace(self):
        return self.grad().divergence()
    def simplify(self):
        return self
    def __repr__(self):
        return (f"General object '{self.name}', attached objects: {self.attached}")
    def __str__(self):
        return f"{self.name}({self.attached})"
    def __eq__(self,other):
        if not isinstance(other,GeneralObject):
            return False
        else:
            if isinstance(other,Complex) and other.isreal():
                return self==other.real
            else:
                return self.name==other.name and self.attached==other.attached

    def __getitem__(self,sl):
        if isinstance(sl,int):
            try:
                return self.attached[sl-1]
            except IndexError: 
                if isinstance(self,Vector):
                    raise ValueError(f"The vector doesn't have a component with index {sl}!")
                elif isinstance(self,Complex):
                    raise ValueError(f"The complex number doesn't have a component with index {sl}!")
                else:
                    raise ValueError(f"The expression {self} doesn't have an argument with index {sl}!")
        if isinstance(sl,Vector):
            try: 
                if all(isinstance(comp,Variable) for comp in sl):
                    return ft.reduce(lambda x,y: x.deriv(y),sl,self)
                #print("vec:",sl)
                sl=slice(*sl)
                #print("slice:",sl)
                return Vector(*self.attached[slice(sl.start-1 if not (sl.start is  None or sl.start==0) else None,sl.stop  if not (sl.stop is None or sl.stop==0) else None,sl.step)])
            except (AttributeError,TypeError) as e:
                raise ValueError(f"The indices must either all be integers or all be variables ! |{e}")
        elif isinstance(sl,Variable):
            return self.deriv(sl)
        elif isinstance(sl,slice):
            return Vector(*self.attached[slice(sl.start-1 if not (sl.start is  None or sl.start==0) else None,sl.stop  if not (sl.stop is None or sl.stop==0) else None,sl.step)])
        else:
            raise ValueError("The index must either be an integer or a variable or a vector of integers or a vector of variables")
    
    def __iter__(self):
        raise TypeError(f"The object {self.__class__} is not iterable")                                                                                                 
    
    def __hash__(self):
        return hash(repr(self))
    
    @property
    def numpyfunc(self):
        return NotImplemented
        """def undef_func(*args,**kwargs):
            raise NotImplementedError(f"Type {self.__class__} doesn't implement numpyfunc!")
        
        return undef_func
        """

    def numpiable(self):
        if self.numpyfunc is NotImplemented:
            return False
        else:
            return all(not isinstance(arg,GeneralObject) or arg.numpiable() for arg in self.attached)
    
    @property
    def numcompose(self):
        if not self.numpiable():
            raise TypeError("This object can't be numpified!")
        if isinstance(self,Variable):
            return self.numpyfunc
        funcgen=[a.numpyfunc if isinstance(a,Variable) 
                else a.numcompose if isinstance(a,GeneralObject) 
                else numpyconst(a)
                for a in self.attached]
        def res(**kwargs):
            return self.numpyfunc(*(func(**kwargs) for func in funcgen))
        return res

class Variable(GeneralObject):
    @property
    def name(self):
        return "Variable"
    @property 
    def evaluables(self):
        return (int,float,complex)
    @property
    def is_evaluable(self):
        return all(obj.is_evaluable for obj in self.attached)
    def __init__(self,varname):
        self.attached=varname,
        self.varname=varname

    def eval(self,*args,**kwargs):
        for name,value in kwargs.items():
            if name==self.varname:
                return value
        return self

    def numpyfunc(self,*args,**kwargs):
        #print(kwargs)
        res=self.eval(*args,**kwargs)
        if not isinstance(res,np.ndarray):
            for var,val in kwargs.items():
                if isinstance(val,np.ndarray):
                    res=res*np.ones(np.shape(val))
                    break
        #print(kwargs,res)
        
        if isinstance(res,GeneralObject):
            
            if isinstance(res,Complex):
                try:
                    return complex(res)
                except TypeError:
                    raise ValueError("The expression can't contain any indeterminate variables!")
            elif isinstance(res,Vector):
                res=list(res.flat())
                rl=[]
                for comp in res:
                    if isinstance(comp,GeneralObject):
                        if isinstance(comp,Complex):
                            try:
                                rl.append(complex(comp))
                            except TypeError:
                                raise ValueError("None of the components of the final vector may contain any indeterminate variables!")
                        else:
                            raise ValueError("None of the components of the final vector may contain any indeterminate variables!")
                    else:
                        rl.append(comp)
                
                return np.array(rl)
            else:
                raise ValueError("The expression can't contain any indeterminate variables!")
        else:
            return res
        
    
    def deriv(self,var):
        if var==self.varname or var==self:
            return Monom(1) #1
        else:
            return Monom(0) #0
    def __repr__(self):
        return self.varname+" [Variable]"
    def __str__(self):
        
        return self.varname
    def get_vars(self):
        #print("vr",self.varname)
        return {self.varname}
    def get_const_vars(self):
        #print("vr",self.varname)
        return set()

class PrettyVariable(Variable):
    def __init__(self,varname):
        super().__init__(self.pretty(varname))
    @staticmethod
    def pretty(varname):
        if varname:
            try:
                
                if varname[0].islower():
                    return unicodedata.lookup(f"Greek small letter {varname}")
                if varname[0].isupper():
                    return unicodedata.lookup(f"Greek Capital letter {varname}")
                else:
                    return varname 
            except KeyError:
                try:
                    return unicodedata.lookup(f"Hebrew letter {varname}")
                except KeyError:
                    return varname
        else:
            raise ValueError("The variable must have a name!")

class Const(GeneralObject):
    @property
    def name(self):
        return "const"
    @property 
    def evaluables(self):
        return (int,float,complex)
    @property
    def is_evaluable(self):
        return all(obj.is_evaluable for obj in self.attached)
    def __init__(self,varname,value):
        self.attached=[varname,value]
        self.varname=varname
        self.value=value
    def eval(self,*args,**kwargs):
        for name,_ in kwargs.items():
            #print(repr(name),repr(self.varname),name==self.varname)
            if name==str(self.varname):
                return self.value
        return self
    @property
    def numpyfunc(self,*args,**kwargs):
        #return NotImplemented
        def func(x,y): return numpyconst(y)(whatever=0)#(*args)*self.value
        #return lambda x,y: numpyconst(y)(asdasad=0)#func
        return func
        #return self.eval(*args,**kwargs)
        
    def deriv(self,var):
        return Monom(0) #0
    def contains_const(self):
        return True        
    def get_const_vars(self):
        return {str(self.varname)}
    def __repr__(self):
        return f"const({self.varname},{self.value})"

    def __str__(self):
        return repr(self)


class VariableFunction(GeneralObject):

    @property
    def name(self):
        return "VariableFunction"
    @property 
    def evaluables(self):
        return (int,float,complex,GeneralObject)
    @property
    def is_evaluable(self):
        return all(obj.is_evaluable for obj in self.attached)

    def __init__(self,varname,*arguments,derorder=None):
        self.funcname=varname
        self.attached=self.arguments=arguments
        if derorder is None:
            self.derorder=[0 for _ in self.arguments]
        else:
           # print(derorder,self.arguments)
            if len(derorder)==len(self.arguments):
               self.derorder=derorder
            else:
                raise ValueError("The dimension of derorder must be equal to the number of functional arguments")
            
        
    def eval(self,*args,**kwargs):
        if self.funcname in kwargs:
            name,value=self.funcname,kwargs[self.funcname]
            if all(order==0 for order in self.derorder):
                return value(*self.arguments).eval(*args,**kwargs)
            else:
                dummyargs=[Variable(f"dummy_{i}") for i,_ in enumerate(self.arguments)]
                func=value(*dummyargs) #self.arguments
                try:
                    for pos,order in enumerate(self.derorder):
                        for _ in range(order):
                            func=func.deriv(dummyargs[pos])
                    func=func.eval(**{darg.varname:arg for darg,arg in zip(dummyargs,self.arguments)})
                    return func.eval(*args,**kwargs)
                except AttributeError:
                    raise NotImplementedError(f"The General object of type {type(func)}  does not specify their derivatives yet!")
                
        else:
            return VariableFunction(self.funcname,*(arg.eval(**kwargs) if isinstance(arg,GeneralObject) else arg for arg in self.arguments),derorder=self.derorder)

    def deriv(self,var):
        return sum(VariableFunction(self.funcname,*self.arguments,derorder=[order+1 if i==ind else order for i,order in enumerate(self.derorder)])*argument.deriv(var)
        for ind,argument in enumerate(self.arguments))
    def partial(self,*nums):
        #print(self.derorder)
        nderorder=list(self.derorder)
        for num in nums:
            if num in [1,len(self.arguments)]:
                nderorder[num-1]+=1
        return VariableFunction(self.funcname,*self.arguments,derorder=nderorder)
    def __repr__(self):
        return f"VariableFunction {self.funcname}({','.join(str(arg) for arg in self.arguments)}), order of derivative = {self.derorder}"
    def __str__(self):
        return f"{self.funcname}({','.join(str(arg) for arg in self.arguments)})"+("" if all(order==0 for order in self.derorder) else f"_{self.derorder}")

class BinOp(GeneralObject):
    @property
    def name(self):
        return "Bin. op"
    @property
    def evaluables(self):
        #return ( (int,float),(pl.Polynomial,int),(pl.Polynomial,float))
        return ( int,float,Complex,Monom,Polynomial,Quotient,Vector)
    @property
    def is_evaluable(self):
        return all(obj.is_evaluable for obj in self.attached)
    @abstractmethod
    def action(self,*args):
        pass
    @abstractproperty
    def symbol(self):
        return ""
    
    def __init__(self,*attached):
        super().__init__(*attached)
        self.obj1=attached[0]
        self.obj2=attached[1]

    def eval(self,*args,**kwargs): 
        obj1,obj2=(obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in self.attached)
        #print("Multiplying...") 
        if isinstance(obj1,self.evaluables) and isinstance(obj2,self.evaluables):
            return self.action(obj1,obj2)
        else:
            return self.__class__(obj1,obj2)
    def __str__(self):
        return self.symbol.join("("+str(el)+")" for el in self.attached)

class ObjSum(BinOp):
    @property
    def name(self):
        return "ObjSum"
    def action(self,*args):
        #print(args[0]+args[1])
        return sum(args)
    @property
    def symbol(self):
        return "+"
    def deriv(self,dervar):
        xder,yder=(obj.deriv(dervar) if isinstance(obj,GeneralObject) else Monom(0) for obj in self.attached)
        return xder+yder
    @property
    def numpyfunc(self):
        def meval(x,y):
                return x+y
        return meval
class ObjMul(BinOp):
    @property
    def name(self):
        return "ObjMul"
    def action(self,*args):
        
        return ft.reduce(mul,args)
    @property
    def symbol(self):
        return "*"
    def deriv(self,dervar):
        x,y=self.attached
        xder,yder=(obj.deriv(dervar) if isinstance(obj,GeneralObject) else Monom(0) for obj in (x,y))
        return xder*y+x*yder

class ObjPow(BinOp):
    @property
    def name(self):
        return "ObjPow"
    def action(self,*args):
        res=pow(args[0],args[1])
        if isinstance(res,complex):
            return Complex(res)
        else:
            return res
    @property
    def symbol(self):
        return "^"
    @property 
    def numpyfunc(self):
        return lambda x,y: x**y
    def eval(self,*args,**kwargs): 
        obj1,obj2=(obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in self.attached)
        #print("Multiplying...") fca
        if isinstance(obj2,self.evaluables):
            if isinstance(obj1,self.evaluables):    
                return self.action(obj1,obj2)
            elif isinstance(obj2,int) and obj2>0:
                return Monom(obj1)**obj2
            else:
                return ObjPow(obj1,obj2)
        else:
            return ObjPow(obj1,obj2)

    def deriv(self,dervar):
        x,y=self.attached
        xder,yder=(obj.deriv(dervar) if isinstance(obj,GeneralObject) else Monom(0) for obj in (x,y))
        return Ln(x)*ObjPow(x,y)*yder+y*ObjPow(x,y-1)*xder
    def __str__(self):
        lp="("+str(self.attached[0])+")" if len(str(self.attached[0]))>1 else str(self.attached[0])
        up="("+str(self.attached[1])+")" #if len(str(self.attached[1]))>1 else str(self.attached[1])
        
        return lp+"^"+up

class Equation(BinOp):

    def __init__(self,lhs,rhs,solver):
        super().__init__(lhs,rhs)
        self.lhs=lhs
        self.rhs=rhs
        self.solver=solver
    @property
    def name(self):
        return "equation"
    def action(self,*args):
        return self.solver # solve
    @property
    def symbol(self):
        return " = "
    def eval(self,*args,format=False,fractions=False,**kwargs):
        lhs,rhs=(obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in self.attached)
        try:
            return self.solver(lhs,rhs,format=format,fractions=fractions)
        except NotImplementedError: #
            #raise
            return self.__class__(lhs,rhs,solver=self.solver)
    
    def __add__(self,other):
        
        if isinstance(other,Equation):
            return self.__class__(self.lhs+other.lhs,self.rhs+other.rhs,solver=self.solver)
        else:
            return self.__class__(self.lhs+other,self.rhs+other,solver=self.solver)
    __radd__=__add__
    def __neg__(self):
        return self.__class__(-self.lhs,-self.rhs,self.solver)
    def __sub__(self,other):
        if isinstance(other,Equation):
            return self.__class__(self.lhs-other.lhs,self.rhs-other.rhs,solver=self.solver)
        else:
            return self.__class__(self.lhs-other,self.rhs-other,solver=self.solver)
    def __rsub__(self,other):
        return -self+other
    
    def __mul__(self,other):
        if isinstance(other,Equation):
            return self.__class__(self.lhs*other.lhs,self.rhs*other.rhs,solver=self.solver)
        else:
            return self.__class__(self.lhs*other,self.rhs*other,solver=self.solver)
    __rmul__=__mul__
    def __truediv__(self,other):
        if isinstance(other,Equation):
            return self.__class__(self.lhs/other.lhs,self.rhs/other.rhs,solver=self.solver)
        else:
            return self.__class__(self.lhs/other,self.rhs/other,solver=self.solver)
    def __rtruediv__(self,other):
        return self.__class__(other/self.lhs,other/self.rhs,self.solver)
    def __pow__(self,other):
        if isinstance(other,Equation):
            return self.__class__(self.lhs**other.lhs,self.rhs**other.rhs,solver=self.solver)
        else:
            return self.__class__(self.lhs**other,self.rhs**other,solver=self.solver)
    __rpow__=__pow__
    def __str__(self):
        return self.symbol.join(str(el) for el in self.attached)
    def simplify(self):
        return Equation(self.lhs.simplify() if isinstance(self.lhs,GeneralObject) else self.lhs,self.lhs.simplify() if isinstance(self.rhs,GeneralObject) else self.rhs,self.solver)

class Inequality(Equation):
    def __init__(self,lhs,rhs,solver,symbol="<"):
            solver=ft.partial(solver,symbol=symbol)
            super().__init__(lhs,rhs,solver)
            self.ineq_symbol=symbol
    @property
    def name(self):
        return "inequality"
    def __str__(self):
        return self.ineq_symbol.join(str(el) for el in self.attached)

class SystemOfEqs(GeneralObject):
    
    def __init__(self,*attached,system_solver):
        #super().__init__(*attached)

        self.eqlist=[]
        self.params=[]
        
        for arg in attached:
            if isinstance(arg,Equation):
                self.eqlist.append(arg)
            elif isinstance(arg,SystemOfEqs):
                for eq in arg:
                    self.eqlist.append(eq)
            else:
                self.params.append(arg)
                #raise ValueError("The expressions between ';' must must be (systems of) equations!")
                
        self.system_solver=system_solver
        self.attached=self.eqlist
    @property
    def name(self):
        return "system of equations"
    def action(self,*args):
        return self.system_solver # solve
    @property
    def evaluables(self):
        return (int,float,GeneralObject)
    @property
    def is_evaluable(self):
        return all(obj.is_evaluable for obj in self.attached)
    @property
    def symbol(self):
        return " ; "
    def eval(self,*args,**kwargs):
        lhs_list,rhs_list=[],[]
        for eq in self.eqlist:
            lhs,rhs=(obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in eq.attached)
            lhs_list.append(lhs)
            rhs_list.append(rhs)
        try:
            
            return self.system_solver(lhs_list,rhs_list)    
        except NotImplementedError: #
            return self.__class__(*(Equation(lhs,rhs,solver=eq.solver) for lhs,rhs,eq in zip(lhs_list,rhs_list,self.eqlist)),system_solver=self.system_solver)
    
    def __str__(self):
        return self.symbol.join(str(el) for el in self.eqlist)
    def __iter__(self):
        return iter(self.eqlist)
    def simplify(self):
        return SystemOfEqs(*(eq.simplify() for eq in self.eqlist),*(par.simplify() if isinstance(par,GeneralObject) else par for par in self.params),self.system_solver)
class Solution(GeneralObject):
    @property
    def name(self):
        return "solution"
    @property
    def evaluables(self):
        return None
    @property
    def is_evaluable(self):
        return False
    def eval(self,*ars,**kwars):
        return Vector(*(sol.eval(*ars,**kwars) for sol in self.solpols))
    def __init__(self,solution,kernel,solvable,variables):
        self.solution=solution
        self.kernel=kernel
        self.solvable=solvable
        self.variables=[var[:-11] if var.endswith("[Variable]") else var for var in variables]
        self.attached=[Vector(*solution),Vector(*(Vector(*k) for k in (kernel))),self.solvable]
        self.solpols,self.kervars=self.make_solpols()
    
    def make_solpols(self):
        sols=[]
        kervars=[]
        for _ in self.kernel:
            for num in list(range(112,122))+list(range(97,122)):
                symb=chr(num)
                if symb not in self.variables+[str(k) for k in kervars]:
                    kervars.append(Variable(symb))
                    break
        for ind,(_,sol) in enumerate(zip(self.variables,self.solution)):
            sol=Monom(float(sol))
            for kvar,k in zip(kervars,self.kernel[:,ind]):
                #print(k,kvar)
                sol+=k*kvar
            sols.append(sol)
        return Vector(*sols),kervars
    def __iter__(self):
        return iter(self.solpols)
    def __repr__(self):
        s=f"Solution object(Solution: {self.solution}, Kernel: {self.kernel}, Exact solution: {self.solvable})"
        return s
    def __str__(self):
        if self.solvable:
            return " , ".join(str(sol) for sol in self.solpols)
        else:
            return "No solution exists, showing L2 approximation: "+" , ".join(str(sol) for sol in self.solpols)
        kervars=[str(kvar) for kvar in self.kervars]
        for _ in self.kernel:
            for num in list(range(115,122))+list(range(97,115)):
                symb=chr(num)
                if symb not in self.variables:
                    kervars.append(symb)
                    break
        def normstr(num):
            if num>0:
                if num==1:
                    return "+"
                else:
                    return "+"+truncate_number(abs(num))
            elif num<0:
                if num==-1:
                    return "-"
                else:
                    return "-"+truncate_number(abs(num))
            else:
                return "+0"
        sparts=[]
        for ind,(var,sol) in enumerate(zip(self.variables,self.solution)):
            s=""
            s+=f"{var} = {truncate_number(sol)}"
            for kvar,k in zip(kervars,self.kernel[:,ind]):
                s+=f"{normstr(k)}{kvar}"
            sparts.append(s)
        return " , ".join(sparts)
    
class Func(GeneralObject):
    
    
    @property
    def evaluables(self):
        return int,float
    @property
    def is_evaluable(self):
        return all(isinstance(obj,self.evaluables) for obj in self.attached )
    @abstractproperty
    def numarg(self):
        return None

    @abstractproperty
    def func(self):
        return (lambda : None)
    @abstractproperty
    def inverse(self):
        return type(None)

    def eval(self,*args,**kwargs): #nefunguje správně, co když některé argumenty budou čísla a jiná ne? Aha, možná to funguje správně
        #print(self.attached)
        #print("report:",str(self),self.attached)
        if len(self.attached)!=self.numarg:
            raise ValueError(f"Func {self.name}  needs to take exactly {self.numarg} argument(s)!")
        else:
            
            #if Vector not in self.evaluables and  isinstance(self.attached[0],Vector):
            #    return Vector(*(self.__class__(comp).eval(*args,**kwargs) for comp in self.attached[0]))
            results=[obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in self.attached]
            #print(self,results,all(isinstance(res,self.evaluables) for res in results))
            if all(isinstance(res,self.evaluables) for res in results):
                return self.func(*(results))
            else:
                if Vector not in self.evaluables:
                    for argind,arg in enumerate(results):
                        if isinstance(arg,Vector):
                            arglist=[ results[:argind]+[comp]+results[argind+1:] for comp in arg]
                            return arg.__class__(*(self.__class__(*arguments).eval(*args,**kwargs) for arguments in arglist))
            #    tr=self.__class__(*results)
                return self.__class__(*results)
    @classmethod
    def complexification(cls,realfunc,complfunc):
        def compeval(z,*args):
            if isinstance(z,(Complex,complex)):
                if isinstance(z.real,(int,float)) and isinstance(z.imag,(int,float)):
                    res=complfunc(z.real+1j*z.imag,*args)
                    return Complex(res.real,res.imag)
                else:
                    return cls(z,*args)
            else:    
                return realfunc(z,*args)
        return compeval

    def __str__(self):
        s=",".join(str(arg) for arg in self.attached)
        return f"{self.name}({s})"

class Exp(Func):
    @property
    def name(self):
        return "exp"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Ln
    @property
    def func(self):
        return self.complexification(math.exp,cmath.exp)
    @property
    def numpyfunc(self):
        return np.exp
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,Polynomial):
            return product((Exp(mon) if mon.coef>=0 else 1/Exp(-mon)  for mon in arg.monlist))
        elif isinstance(arg,Complex):
            if arg.isimag():
                return Complex(Cos(arg.imag),Sin(arg.imag))
            else:
                return (Exp(arg.real)*Complex(Cos(arg.imag),Sin(arg.imag)))
        return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return self*(self.attached[0]).deriv(dervar)

class Ln(Func):
    @property
    def name(self):
        return "ln"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Exp
    @property
    def func(self):
        f=self.complexification(math.log,cmath.log)
        def corlog(z):
            try:
                return f(z)
            except ValueError:
                raise ZeroDivisionError(
            "Can't take logarithm of a non-positive number. For non-zero numbers, change to complex arguments (x->x+0i).")
        return corlog
    @property
    def numpyfunc(self):
        return np.log
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):
            
            return arg.attached[0]
        if isinstance(arg,Monom):
            l1=sum(order*self.__class__(obj).simplify()  for var,(obj,order) in arg.vars.items())
            return l1+self.__class__(arg.coef) if arg.coef!=1 else l1
        elif isinstance(arg,Polynomial) and len(arg.monlist)==1:
            l1= sum(order*self.__class__(obj).simplify()  for var,(obj,order) in arg.monlist[0].vars.items())
            return l1+self.__class__(arg.coef) if arg.coef!=1 else l1
        elif isinstance(arg,ObjPow):
            return arg.attached[1]*self.__class__(arg.attached[0]).simplify()
        elif isinstance(arg,Complex):
            if arg.isimag():
                try:
                    if arg.imag>0:
                        return Complex(math.log(arg.imag),math.pi/2)
                    else:
                        return Complex(math.log(-arg.imag),-math.pi/2)
                except (TypeError,ValueError):
                    return Complex(Ln(arg.norm()),arg.angle()).eval()
            else:
                return Complex(Ln(arg.norm()),arg.angle()).eval()
        else:
            return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/self.attached[0]*(self.attached[0]).deriv(dervar) 


class Log(Func):
    @property
    def name(self):
        return "log"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        f=self.complexification(math.log10,cmath.log10)
        def corlog(z):
            try:
                return f(z)
            except ValueError:
                raise ZeroDivisionError(
            "Can't define take logarithm of a non-positive number. For non-zero numbers, change to complex arguments (x->x+0i). ")
        return corlog
    @property
    def numpyfunc(self):
        return np.log10
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):
            
            return arg.attached[0]
        if isinstance(arg,Monom):
            l1=sum(order*self.__class__(obj).simplify()  for var,(obj,order) in arg.vars.items())
            return l1+self.__class__(arg.coef) if arg.coef!=1 else l1
        elif isinstance(arg,Polynomial) and len(arg.monlist)==1:
            l1= sum(order*self.__class__(obj).simplify()  for var,(obj,order) in arg.monlist[0].vars.items())
            return l1+self.__class__(arg.coef) if arg.coef!=1 else l1
        elif isinstance(arg,ObjPow):
            return arg.attached[1]*self.__class__(arg.attached[0]).simplify()
        
        elif isinstance(arg,Complex):
            if arg.isimag():
                try:
                    if arg.imag>0:
                        return Complex(math.log10(arg.imag),math.pi/2)
                    else:
                        return Complex(math.log10(-arg.imag),-math.pi/2)
                except (TypeError,ValueError):
                    return Complex(Log(arg.norm()),arg.angle()).eval()
            else:
                return Complex(Log(arg.norm()),arg.angle()).eval()
        else:
            return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/self.attached[0]*math.log10(math.e)*(self.attached[0]).deriv(dervar) 
#New code ↓
class Sinh(Func):
    @property
    def name(self):
        return "sinh"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Asinh
    @property
    def func(self):
        return self.complexification(math.sinh,cmath.sinh)
    @property
    def numpyfunc(self):
        return np.sinh
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
                return 2*(Sinh(arg/2)*Cosh(arg/2))
        elif isinstance(arg,Complex):
            if arg.isimag():
                return Complex(0,Sin(arg.imag)).eval()
            else:
                return Complex(Sinh(arg.real)*Cos(arg.imag),Cosh(arg.real)*Sin(arg.imag)).eval()
        else:
            return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Cosh(self.attached[0])*(self.attached[0]).deriv(dervar)

class Asinh(Func):
    @property
    def name(self):
        return "asinh"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Sinh
    @property
    def func(self):
            return self.complexification(math.asinh,cmath.asinh)
    @property
    def numpyfunc(self):
        return np.arcsinh
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,Complex):
            return  Ln(arg+Sqrt(arg**2+1)).eval()
        else:
            return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/Sqrt(1+self.attached[0]**2)*self.attached[0].deriv(dervar)
class Cosh(Func):
    @property
    def name(self):
        return "cosh"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Acosh
    @property
    def func(self):
        return self.complexification(math.cosh,cmath.cosh)
    @property
    def numpyfunc(self):
        return np.cosh
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        
        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
                return Cosh(arg/2)**2+Sinh(arg/2)**2
        elif isinstance(arg,Complex):
            if arg.isimag():
                return Complex(Cos(arg.imag),0).eval()
            else:
                return Complex(Cosh(arg.real)*Cos(arg.imag),Sinh(arg.real)*Sin(arg.imag)).eval()
        else:
            return self

    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Sinh(self.attached[0])*(self.attached[0]).deriv(dervar)
class Acosh(Func):
    @property
    def name(self):
        return "acosh"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Cosh
    
    @property
    def func(self):
        return self.complexification(math.acosh,cmath.acosh)
    @property
    def numpyfunc(self):
        return np.arccosh
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,Complex):
            return Ln(arg+Sqrt(arg**2-1)).eval()
        return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/Sqrt(-1+self.attached[0]**2)*self.attached[0].deriv(dervar)
#New code ↑
class Gamma(Func):
    @property
    def name(self):
        return "Gamma"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        ggamma=self.complexification(math.gamma,cgamma)
        def ext_func(z):
            try:
                res=ggamma(z)
                try:
                    ires=int(res)
                except TypeError:
                    return res
                return ires if ires==res else res
            except ValueError:
                raise ZeroDivisionError("The gamma function is undefined for negative integers!")
          #  try:
          #      return math.factorial(x)
          #  except ValueError:
          #      return ggamma(x+1)
        return ext_func
    @property
    def numpyfunc(self):
        return cgamma
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):            
            return arg.attached[0]

        return self
    def deriv(self,dervar):
        raise NotImplementedError("The derivative of factorial/Gamma function is not implemented yet")
    def __str__(self):
        return "("+str(self.attached[0]-1)+")!"

class Beta(Func):
    @property
    def name(self):
        return "betafunc"
    @property
    def numarg(self):
        return 2
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def cbeta(x,y):
            if any(isinstance(a,Complex) for a in (x,y)):
                return (Gamma(x)*Gamma(y)/(Gamma(x+y))).eval()    
            else:
                return beta(x,y)
        return cbeta
      
    @property
    def numpyfunc(self):
        def npbeta(x,y):
            return cgamma(x)*cgamma(y)/cgamma(x+y)
        return npbeta
    def simplify(self):
        return self
    def deriv(self,dervar):
        raise NotImplementedError("The derivative of factorial/Gamma function is not implemented yet")

    
class Erf(Func):
    @property
    def name(self):
        return "erf"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        return self.complexification(erf,erf)
    @property
    def numpyfunc(self):
        return erf
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):            
            return arg.attached[0]

        return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 2/math.sqrt(math.pi)*Exp(-self.attached[0]**2)*self.attached[0].deriv(dervar)
 #   def __str__(self):
  #      return "("+str(self.attached[0]-1)+")!"

class Erfc(Func):
    @property
    def name(self):
        return "erfc"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        return self.complexification(erfc,erfc)
    @property
    def numpyfunc(self):
        return erfc
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):            
            return arg.attached[0]

        return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return -2/math.sqrt(math.pi)*Exp(-self.attached[0]**2)*self.attached[0].deriv(dervar)
    #def __str__(self):
    #    return "("+str(self.attached[0]-1)+")!"

class OldAbs(Func):
    @property
    def name(self):
        return "abs"
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        return abs
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):
            
            return arg.attached[0]
        
        return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return (self.attached[0])/Abs(self.attached[0])*(self.attached[0]).deriv(dervar)
class Sqrt(Func):
    @property
    def name(self):
        return "sqrt"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex,complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        return self.complexification(math.sqrt,cmath.sqrt)
    @property
    def numpyfunc(self):
        return np.sqrt
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):    
            return arg.attached[0]
        if isinstance(arg,Complex):
            return (arg**(1/2)).eval()
        elif isinstance(arg,(int,float)) and arg<0:
            return Complex(0,math.sqrt(-arg))
        else:
            return self

    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0.5*1/Sqrt(self.attached[0])*(self.attached[0]).deriv(dervar)
    def __str__(self):
        return "√("+str(self.attached[0])+")"
class Sin(Func):
    @property
    def name(self):
        return "sin"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Asin

#newcode ↓
    @property
    def func(self):
        return self.complexification(math.sin,cmath.sin)
    @property
    def numpyfunc(self):
        return np.sin
#newcode ↑

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,Acos):
            return Sqrt(1-arg.attached[0]**2)
        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
                return 2*(Sin(arg/2)*Cos(arg/2))
        if isinstance(arg,Complex):
            if arg.isimag():
                return Complex(0,Sinh(arg.imag)).eval()
            else:
                return Complex(Sin(arg.real)*Cosh(arg.imag),Cos(arg.real)*Sinh(arg.imag)).eval()
        else:
            return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        #print("Sin:",type(Cos(self.attached[0])),type(self.attached[0].deriv(dervar)))
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Cos(self.attached[0])*self.attached[0].deriv(dervar)
    @property
    def partials(self):
        return [Cos(*self.attached)]
class Asin(Func):
    @property
    def name(self):
        return "asin"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Sin
    @property
    def func(self):
        return self.complexification(math.asin,cmath.asin)
    @property
    def numpyfunc(self):
        return np.arcsin
    def simplify(self):
        arg=self.attached[0]
        
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/Sqrt(1-self.attached[0]**2)*self.attached[0].deriv(dervar)
   
class Cos(Func):
    @property
    def name(self):
        return "cos"
    
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Acos
#newcode ↓
    @property
    def func(self):
        return self.complexification(math.cos,cmath.cos)
    @property
    def numpyfunc(self):
        return np.cos
#newcode ↑
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,Asin):
            return Sqrt(1-arg.attached[0]**2)
        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
            return Cos(arg/2)**2-Sin(arg/2)**2
        elif isinstance(arg,Complex):
            if arg.isimag():
                return Complex(Cosh(arg.imag),0).eval()
            else:
                return Complex(Cos(arg.real)*Cosh(arg.imag),-Sin(arg.real)*Sinh(arg.imag)).eval()
        else:
            return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return -Sin(*self.attached)*self.attached[0].deriv(dervar)
    @property
    def partials(self):
        return [-Sin(*self.attached)]

class Acos(Func):
    @property
    def name(self):
        return "acos"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Cos
    @property
    def func(self):
        return self.complexification(math.acos,cmath.acos)
    @property
    def numpyfunc(self):
        return np.arccos
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return self
    def deriv(self,dervar):
        #print(self,self.attached[0])
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return -1/Sqrt(1-self.attached[0]**2)*self.attached[0].deriv(dervar)

class Tan(Func):
    @property
    def name(self):
        return "tan"
    
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return Atan
    @property
    def evaluables(self):
        return int,float,Complex
    @property
    def func(self):
        f=self.complexification(math.tan,cmath.tan)
        def cortan(z):
            try:
                return f(z)
            except ValueError:
                raise ZeroDivisionError

        return cortan
    @property
    def numpyfunc(self):
        return np.tan
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
                return 2*Tan(arg/2)/(1-Tan(arg/2)**2)#**2)
        else:
            return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return -1/Cos(arg)**2*arg.deriv(dervar)

class Atan(Func):
    @property
    def name(self):
        return "atan"
    
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return Tan
    @property
    def evaluables(self):
        return int,float,Complex
    @property
    def func(self):
        return self.complexification(math.atan,cmath.atan)
    @property
    def numpyfunc(self):
        return np.arctan
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]

        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
                return 2*Tan(arg/2)/(1-Tan(arg/2))#**2)
        else:
            return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/(arg**2+1)*arg.deriv(dervar)
class Tanh(Func):
    @property
    def name(self):
        return "tanh"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Atanh
    @property
    def func(self):
        return self.complexification(math.tanh,cmath.tanh)
    @property
    def numpyfunc(self):
        return np.tanh
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        
        
        return Sinh(arg)/Cosh(arg)
        #return self

    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return (1-Tanh(self.attached[0])**2)*(self.attached[0]).deriv(dervar)
class Atanh(Func):
    @property
    def name(self):
        return "atanh"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return Tanh
    @property
    def func(self):
        return self.complexification(math.atanh,cmath.atanh)
    @property
    def numpyfunc(self):
        return np.arctanh
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return 0.5*Ln( (1+self.attached[0])/(1-self.attached[0]))
        #return self

    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 1/(1-self.attached[0]**2)*(self.attached[0]).deriv(dervar)
class Real(Func):
    @property
    def name(self):
        return "Real"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,complex,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def tr(val):
            if isinstance(val,(Complex,complex)):
                return val.real
            else:
                return val
        return tr
    @property
    def numpyfunc(self):
        return np.real
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,Complex):
            return arg.real
        return self
        #return self

    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Real(self.attached[0].deriv(dervar))
class Imag(Func):
    @property
    def name(self):
        return "Imag"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,complex,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def tr(val):
            if isinstance(val,(Complex,complex)):
                return val.imag
            else:
                return Monom(0)
        return tr
    @property
    def numpyfunc(self):
        return np.imag
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,Complex):
            return arg.imag
        return self
        #return self

    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Imag(self.attached[0].deriv(dervar))
class Bessel(Func):
    @property 
    def name(self):
        return "J1"
    @property
    def evaluables(self):
        return int,float,complex,Complex
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def fast_bessel(z,n):
            if int(n)==0:
                return scipy.special.j0(z)
            elif int(n)==1:
                return scipy.special.j1(z)
            else:
                return scipy.special.jv(n,z)
        def bessel(z,n):
            return scipy.special.jv(n,z)
        return self.complexification(fast_bessel,bessel)
    @property
    def numpyfunc(self):
        def bessel(z,n):
            return scipy.special.jv(n,z)
        return bessel

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return self

    def deriv(self,dervar):
        z,n=self.attached

        if isinstance(z,(int,float)):
            return 0
        else:
            return 0.5*((self.__class__(n-1,z)-self.__class__(n+1,z)))*z.deriv(dervar)
    def __str__(self):
        s=f"{self.attached[1]},{self.attached[0]}"
        return f"{self.name}({s})"

class SecondBessel(Bessel):
    @property 
    def name(self):
        return "J2"
    @property
    def func(self):
        def fast_bessel(z,n):
            if int(n)==0:
                return scipy.special.y0(z)
            elif int(n)==1:
                return scipy.special.y1(z)
            else:
                return scipy.special.yv(n,z)
        def bessel(z,n):
            return scipy.special.yv(n,z)
        return self.complexification(fast_bessel,bessel)
    @property
    def numpyfunc(self):
        def second_bessel(z,n):
            return scipy.special.yv(n,z)
        return second_bessel

class Ellip1(Func):
    @property 
    def name(self):
        return "E1"
    @property
    def evaluables(self):
        return int,float,complex,Complex
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def ellip1(phi,k):
            return ellipkinc(phi,k**2)
        return ellip1
    @property
    def numpyfunc(self):
        def ellip1(phi,k):
            return ellipkinc(phi,k**2)
        return ellip1

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return self

    def deriv(self,dervar):
        raise NotImplementedError("The derivative of the elliptic integrals are not implemented yet")

class Ellip2(Func):
    @property 
    def name(self):
        return "E2"
    @property
    def evaluables(self):
        return int,float,complex,Complex
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def ellip2(phi,k):
            return ellipeinc(phi,k**2)
        return ellip2
    @property
    def numpyfunc(self):
        def ellip2(phi,k):
            return ellipeinc(phi,k**2)
        return ellip2

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return self

    def deriv(self,dervar):
        raise NotImplementedError("The derivative of the elliptic integrals are not implemented yet")

class CompEllip1(Ellip1):
    @property
    def name(self):
        return "CE1"
    @property
    def numarg(self):
        return 1
    @property
    def func(self):
        def ellip1(phi,k):
            return ellipkinc(np.pi/2,k**2)
        return ellip1
    @property
    def numpyfunc(self):
        def ellip1(k):
            return ellipkinc(np.pi/2,k**2)
        return ellip1

class CompEllip2(Ellip2):
    @property
    def name(self):
        return "CE2"
    @property
    def numarg(self):
        return 1
    @property
    def func(self):
        def ellip2(phi,k):
            return ellipeinc(np.pi/2,k**2)
        return ellip2
    @property
    def numpyfunc(self):
        def ellip2(k):
            return ellipeinc(np.pi/2,k**2)
        return ellip2
        
class Choose(Func):
    @property 
    def name(self):
        return "choose"
    @property
    def evaluables(self):
        return int,float,Complex
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def choose(a1,a2):
            try:
                return (Gamma(a1+1)/(Gamma(a2+1)*Gamma(a1-a2+1))).eval()
            except ZeroDivisionError:
                raise ZeroDivisionError("Can't evalute choose(n,k) for n,k integers and n<k.")
            except OverflowError:
                return binom(a1,a2)
        return choose
    @property
    def numpyfunc(self):
        def choose(x,y):
            return cgamma(x)/(cgamma(y)*cgamma(x-y))
        return choose
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        if isinstance(arg,(Monom,Polynomial)) and arg.coef%2==0:
                return 2*Tan(arg/2)/(1-Tan(arg/2))#**2)
        return self
    def deriv(self,dervar):
        raise NotImplementedError("The derivative of factorial/Gamma function is not implemented yet")

class BayesianFormula(Func):
    @property 
    def name(self):
        return "BT"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def func(proportion,sensitivity,specificity):
            try:
                #if not all(0<=x<=1 for x in (proportion,sensitivity,specificity)):
                #    raise ValueError("All arguments (proportion, sensitivity and specificity) must be real numbers between 0 and 1")
                num=sensitivity*proportion
                return num/(num+(1-proportion)*(1-specificity))
            except ZeroDivisionError:
                if proportion==0:
                    return 0
                else:
                    return 1
        return func
    @property
    def numpyfunc(self):
        def func(proportion,sensitivity,specificity):
            try:
                #if not all( all((0<=x)*(x<=1)) for x in (proportion,sensitivity,specificity)):
                #    raise ValueError("All arguments (proportion, sensitivity and specificity) must be real numbers between 0 and 1")
                num=sensitivity*proportion
                return num/(num+(1-proportion)*(1-specificity))
            except ZeroDivisionError:
                return np.where(proportion==0,0,1)
        return func
    def simplify(self):
        return self
    def deriv(self,dervar):
        x,y,z=self.attached
        denom=x*y+(1-x)*(1-z)
        expr=Quotient(x*y,denom)
        return expr.deriv(dervar)

class RandomSample(Func): #Popřípadě ošetřit další operace:
    @property
    def name(self):
        return "rand"
    @property 
    def evaluables(self):
        return (int,float,complex)
    @property
    def numarg(self):
        return 1
    def get_vars(self):
        return self.vars

    def __init__(self,n=1,x0=0,x1=1,func=1):
        if isinstance(n,int):
            self.is_func=False
            self.n=n
        if isinstance(n,GeneralObject):
            self.is_func=True
            self.vars=n.get_vars()
            self.n=len(self.vars)
            
        if not isinstance(func,GeneralObject):
            
            func=Const(Variable("x"),1)
        if x0==0 and x1==-1:
            self.lims="uni"
        else:
            try:
                iter(x0)
            except TypeError:
                try:
                    iter(x1)
                except TypeError:
                    self.lims=np.array([x0,x1]).reshape(1,2)
            else:
                self.lims=np.array([(x0i,x1i) for x0i,x1i in zip(x0,x1) ])

        self.f=func
        self.x0=x0;self.x1=x1
        self.attached=[n]
        self.samp_vars=func.get_vars()
        self.dim=len(self.samp_vars)
        
    def eval(self,*args,**kwargs):
        if not self.is_func or set(kwargs.keys())==self.vars:
            return super().eval(*args,**kwargs)
        else:
            return self
    @property
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def genrand(x):
            n=5*(1000+self.n)
            #print(self.f,n,self.lims,self.dim)
            res=self.mh(self.f,n=n,lims=self.lims,dim=self.dim)[n//5-1000::4][1000:]
            if self.dim==1:
                res=res[:,0]
            else:
                res=ImplicitVector(*(Vector(*comp) for comp in res[-self.n:]))
            if self.n==1:
                return res[-1]
            else:
                return ImplicitVector(*res[-self.n:])
        return genrand
    @property
    def numpyfunc(self,*args,**kwargs):
        #return NotImplemented
        def npfunc(x):
            sh=x.shape
            n=ft.reduce(mul,sh)
            
            n=5*(1000+n)

            res=self.mh(self.f,n=n,lims=self.lims,dim=self.dim)[n//5::4][1000:]
            if self.dim==1:
                return res[:,0].reshape(*sh)
            else:
                return res.reshape(self.dim,*sh)
        return npfunc
        
    def deriv(self,var):
        return Monom(0) #0
    def __repr__(self):
        return f"{self.name}"+"("+",".join(str(comp) for comp in [self.n,self.x0,self.x1,self.f])+")"
    def __str__(self):
        return repr(self)
    def mh(self,f,n=1,dim=1,lims="uni"):
        unbounded=False
        if lims is "uni":
            lims=np.tile([0,1],(dim,1))
            unbounded=True
        lens=[lim[1]-lim[0] for lim in lims ]
        if type(f)==Const:
            return (lims[:,0]+np.random.random(n*dim)*lens).reshape(n,dim)
        elif type(f)==NormalDist:
            _,μ,σ=f.attached
            res=np.random.normal(μ,σ,n*dim)
            if unbounded:
                return res.reshape(n,dim)
            else:
                res=np.where(res<lims[:,0],lims[:,0],res)
                return np.where(res>lims[:,1],lims[:,1],res).reshape(n,dim)
        x=lims[:,0]+np.random.random(dim)*lens
        samples=[]
        #print(f,self.samp_vars)
        fc=f.eval(**{var:x[i] for i,var in enumerate(self.samp_vars)})
        for _ in range(n):
            xn=lims[:,0]+np.random.random(dim)*lens
            fn=f.eval(**{var:xn[i] for i,var in enumerate(self.samp_vars)})
            
            if fn>=fc:
                x=xn
                fc=fn
            else:
                p=np.random.random()
                a=fn/fc
                if p<a:
                    x=xn
                    fc=fn
            samples.append(x)
        return np.array(samples)

    """
    def __add__(self,x):
        if isinstance(x,Vector):
            return ObjSum(self,x)
        else:
            return super().__add__(x)
    def __radd__(self,x):
        if isinstance(x,Vector):
            return ObjSum(x,self)
        else:
            return super().__radd__(x)
    """
class VectorFunc(Func):
    @property 
    def name(self):
        return "none"
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def evaluables(self):
        return Vector,int,float,complex
    def eval(self,*args,**kwargs): #nefunguje správně, co když některé argumenty budou čísla a jiná ne? Aha, možná to funguje správně
        #print(self.attached)
        if len(self.attached)!=self.numarg:
            raise ValueError(f"Func {self.name}  needs to take exactly {self.numarg} argument(s)!")
        else:
            if any(isinstance(el,(int,float)) for el in self.attached):
                raise ValueError("The arguments must be variable or tensorial expressions")
            if Vector not in self.evaluables and  isinstance(self.attached[0],Vector):
                return Vector(*(self.__class__(comp).eval(*args,**kwargs) for comp in self.attached[0]))
            results=[obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in self.attached]
            if all(isinstance(res,self.evaluables) for res in results):
                return self.func(*(results))
            else:
            #    tr=self.__class__(*results)
                return self.__class__(*results)
    @property
    def func(self):
        return NotImplemented
    @property
    def numpyfunc(self):
        return NotImplemented
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return self
    def deriv(self,dervar):
        v1,v2=self.attached
        return self.__class__(v1.deriv(dervar) if isinstance(v1,GeneralObject) else Monom(0) ,v2)+self.__class__(v1,v2.deriv(dervar) if isinstance(v2,GeneralObject) else Monom(0)) #Tohle bude řádně fungovat až teprve s "obecnou" derivací
    
class Dot_product(VectorFunc):
    @property 
    def name(self):
        return "dot"
    @property
    def numarg(self):
        return 2
    @property
    def func(self):
        def temp_dot(a1,a2): return a1.dot(a2)
        return temp_dot
    @property
    def numpyfunc(self):
        def npdot(xvectors,yvectors):
            xdim,ydim=xvectors.ndim,yvectors.ndim
            if xdim==ydim==1:
                return np.array([vec1.dot(vec2) for vec1,vec2 in zip(xvectors,yvectors)])
            elif xdim==ydim==2:
                sh=xvectors.shape
               # print(xvectors)
                return np.array([vec1.dot(vec2) for xrow,yrow in zip(xvectors,yvectors) for vec1,vec2 in zip(xrow,yrow)  ]).reshape(sh)
            else:
                raise NotImplementedError("Plotting dot product is not implemented for dim>2")
            #for c1,c2 in zp(x,y):
            #    return c1*c2
        return npdot

class FourierTransform(VectorFunc):
    @property 
    def name(self):
        return "FT"
    @property 
    def inverse(self):
        return InverseFourierTransform
    @property
    def numarg(self):
        return 1
    @property
    def func(self):
        def naive_FT(a): 
            N=len(a)
            return 1/N*Vector(*[Complex(sum(a[i+1]*np.exp(-2*np.pi*1j/N*i*k) for i in range(N) )) for k in range(N)])
        def DFT(a):
            n=len(a)
            if (n & (n-1) == 0) and n != 0:
                res=FFT(a)
                return Vector(*(Complex(comp) for comp in res))
            else:
                return naive_FT(a)
        
        return DFT
        #return naive_FT
    @property
    def numpyfunc(self):
        def FFT(a):
            return np.fft.fft(a)
        return FFT
class InverseFourierTransform(VectorFunc):
    @property 
    def name(self):
        return "IFT"
    @property 
    def inverse(self):
        return FourierTransform
    @property
    def numarg(self):
        return 1
    @property
    def func(self):
        def naive_FT(a): 
            N=len(a)
            return Vector(*[Complex(sum(a[i+1]*np.exp(2*np.pi*1j/N*i*k) for i in range(N) )) for k in range(N)])
        def IDFT(a):
            n=len(a)
            if (n & (n-1) == 0) and n != 0:
                res=FFT(a,inverse=True)
                return Vector(*(Complex(comp) for comp in res))
            else:
                return naive_FT(a)
        return IDFT
        #return naive_FT
    @property
    def numpyfunc(self):
        def FFT(a):
            return np.fft.ifft(a)
        return FFT
class Sum(VectorFunc):
    @property 
    def name(self):
        return "sum"
    @property
    def numarg(self):
        return 1
    @property
    def func(self):
        return sum
    @property
    def numpyfunc(self):
        return np.sum #asi bude třeba ještě doplnit osu
class Mul(VectorFunc):
    @property 
    def name(self):
        return "mul"
    @property
    def numarg(self):
        return 1
    @property
    def func(self):
        def mul_func(vec):
            return ft.reduce(mul,vec)
        return mul_func
    @property
    def numpyfunc(self):
        def mul_func(vec):
            return ft.reduce(np.multiply,vec)
        return mul_func #asi bude třeba ještě doplnit osu

class Cross_product(VectorFunc):
    @property 
    def name(self):
        return "cross"
    
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def temp_cross(a1,a2): return a1.cross(a2)
        return temp_cross
class Matmul(VectorFunc):
    @property 
    def name(self):
        return "matmul"
    
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def temp_func(a1,a2): return a1.matmul(a2)
        return temp_func

class Transpose(VectorFunc):
    @property 
    def name(self):
        return "trans"
    
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return self.__class__
    @property
    def func(self):
        def temp_func(a1): return a1.T()
        return temp_func

    def deriv(self,dervar):
        v1=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return self.__class__(v1.deriv(dervar)) #Tohle bude řádně fungovat až teprve s "obecnou" derivací
class Divergence(VectorFunc):
    @property 
    def name(self):
        return "div"
    
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return self.__class__
    @property
    def func(self):
        def temp_func(a1): return a1.divergence()
        return temp_func

    def deriv(self,dervar):
        v1=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return self.__class__(v1.deriv(dervar)) #Tohle bude řádně fungovat až teprve s "obecnou" derivací

class Curl(VectorFunc):
    @property 
    def name(self):
        return "curl"
    
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return self.__class__
    @property
    def func(self):
        def temp_func(a1): return a1.curl()
        return temp_func

    def deriv(self,dervar):
        v1=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return self.__class__(v1.deriv(dervar)) #Tohle bude řádně fungovat až teprve s "obecnou" derivací
class VectorLaplace(VectorFunc):
    @property 
    def name(self):
        return "veclap"
    
    @property
    def numarg(self):
        return 1
    @property 
    def inverse(self):
        return self.__class__
    @property
    def func(self):
        def temp_func(a1): return a1.laplace()
        return temp_func

    def deriv(self,dervar):
        v1=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return self.__class__(v1.deriv(dervar)) #Tohle bude řádně fungovat až teprve s "obecnou" derivací

class Gcd(Func):
    @property 
    def name(self):
        return "gcd"
    
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def evaluables(self):
        return int,
    @property
    def func(self):
        def aux(a1,a2):
            if a1<=0 or a2<=0:
                raise ValueError("All numbers must be positive integers!")
            return gcd(a1,a2)
            
        return aux
    @property
    def numpyfunc(self):
        return NotImplemented
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0
class Lcm(Func):
    @property 
    def name(self):
        return "lcm"
    
    @property
    def numarg(self):
        return 2
    @property 
    def inverse(self):
        return type(None)
    @property
    def evaluables(self):
        return int,
    @property
    def func(self):
        def aux(a1,a2):
            if a1<=0 or a2<=0:
                raise ValueError("All numbers must be positive integers!")
            return lcm(a1,a2)
            
        return aux
    @property
    def numpyfunc(self):
        return NotImplemented
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0
class Sign(Func):
    @property 
    def name(self):
        return "sign"
    
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex,complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def sgn(num):
            try:
                return 1 if num>0 else (-1 if num<0 else 0)
            except TypeError:
                num=Complex(num)
                if num.isreal():
                    return sgn(num.real)
                else:
                    raise ValueError("Can' take sign of a complex number")
        return sgn
    @property
    def numpyfunc(self):
        return np.sign
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0*arg.deriv(dervar) #1/(arg**2+1)*arg.deriv(dervar) Ve skutecnosti by tu mela byt delta funkce  - implementovat?

class NormalDist(Func):
    @property 
    def name(self):
        return "normal"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        c=math.sqrt(2*math.pi)
        def f(x,μ,σ):
            return 1/(σ*c)*Exp(-0.5*((x-μ)/σ)**2).eval()
        return f

    @property
    def numpyfunc(self):
        c=math.sqrt(2*math.pi)
        def f(x,μ,σ):
            return 1/(σ*c)*np.exp(-0.5*((x-μ)/σ)**2)
        return f

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        x,μ,σ=self.attached
        pref=(x-μ)/σ
        return NormalDist(x,μ,σ)*(-pref/σ*(x.deriv(dervar) if isinstance(x,GeneralObject) else 0)+
                                  pref/σ*(μ.deriv(dervar) if isinstance(μ,GeneralObject) else 0)+
                                  1/σ*(pref**2-1)*(σ.deriv(dervar) if isinstance(σ,GeneralObject) else 0))

class LogNormal(Func):
    @property 
    def name(self):
        return "lognormal"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        c=math.sqrt(2*math.pi)
        def f(x,μ,σ):
            return 1/(σ*c*x)*Exp(-0.5*((Ln(x).eval()-μ)/σ)**2).eval()
        return f

    @property
    def numpyfunc(self):
        c=math.sqrt(2*math.pi)
        def f(x,μ,σ):
            return 1/(σ*c*x)*np.exp(-0.5*((np.log(x)-μ)/σ)**2)
        return f

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        x,μ,σ=self.attached
        pref=(Ln(x)-μ)/σ
        return LogNormal(x,μ,σ)*(1/x*(-1-pref)*(x.deriv(dervar) if isinstance(x,GeneralObject) else 0)+
                                  pref*(μ.deriv(dervar) if isinstance(μ,GeneralObject) else 0)+
                                  1/σ*(pref**2-1)*(σ.deriv(dervar) if isinstance(σ,GeneralObject) else 0))

class Maxwell_Boltzmann(Func):
    @property 
    def name(self):
        return "MB"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        
        def f(v,T,M):
            σ=np.sqrt(1000*8.31446261815324*T/M)
            return 2/σ**2*v**2*NormalDist(v,0,σ).eval() 
        return f

    @property
    def numpyfunc(self):
        def f(v,T,M):
            σ=np.sqrt(1000*8.31446261815324*T/M)
            return np.sqrt(2/np.pi)/(σ**3)*v**2*np.exp(-0.5*(v/σ)**2)
        return f

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        v,T,M=self.attached
        σ=Sqrt(1000*8.31446261815324*T/M)
        return (2/σ**2*v**2*NormalDist(v,0,σ)).deriv(dervar) 
                                  
class Binomial(Func):
    @property 
    def name(self):
        return "binom"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        
        def f(p,n,k):

            if n>800 and abs(n-int(n))<self.tol and abs(n-int(n))<self.tol:
                return large_binom(p,int(n),int(k))
            else:
                return (Choose(n,k)*p**k*(1-p)**(n-k)).eval()

        return f

    @property
    def numpyfunc(self):
        def f(p,n,k):
            return binom(n,k)*p**k*(1-p)**(n-k)
        return f

    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        p,n,k=self.attached
        return Monom(0) #LogNormal(x,μ,σ)*(1/x*(-1-pref)*(x.deriv(dervar) if isinstance(x,GeneralObject) else 0)+
                         #         pref*(μ.deriv(dervar) if isinstance(μ,GeneralObject) else 0)+
                          #        1/σ*(pref**2-1)*(σ.deriv(dervar) if isinstance(σ,GeneralObject) else 0))
class BetaDist(Func):
    @property 
    def name(self):
        return "betadist"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def f(x,α,β):
            try:
                return (x**(α-1)*(1-x)**(β-1)/Beta(α,β)).eval()
            except ValueError:
                raise ZeroDivisionError("Can't evaluate the gamma distribution for negative integer α")
        return f

    @property
    def numpyfunc(self):
        def f(x,α,β):
            return x**(α-1)*(1-x)**(β-1)*cgamma(α+β)/(cgamma(α)*cgamma(β))
        return f

    def simplify(self):
        return self
    def deriv(self,dervar):
        raise NotImplementedError("The derivative of the Gamma distribution is not implemented yet")
class GammaDist(Func):
    @property 
    def name(self):
        return "gammadist"
    
    @property
    def numarg(self):
        return 3
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        ggamma=self.complexification(math.gamma,cgamma)
        def f(x,α,β):
            try:
                return β**α*x**(α-1)*np.exp(-β*x)/ggamma(α)
            except ValueError:
                raise ZeroDivisionError("Can't evaluate the gamma distribution for negative integer α")
        return f

    @property
    def numpyfunc(self):
        def f(x,α,β):
            return β**α*x**(α-1)*np.exp(-β*x)/cgamma(α)
        return f

    def simplify(self):
        return self
    def deriv(self,dervar):
        raise NotImplementedError("The derivative of the Gamma distribution is not implemented yet")

class Poisson(Func):
    @property
    def name(self):
        return "poisson"
    @property
    def numarg(self):
        return 2
    @property
    def evaluables(self):
        return int,float,Complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def poisson(k,λ):
            return (λ**k*Exp(-λ)/Gamma(k+1)).eval()
            
        return poisson
    @property
    def numpyfunc(self):
        def nppoisson(k,λ):
            return (λ**k*np.exp(-λ)/cgamma(k+1))
        return nppoisson
    def simplify(self):
        
        return self
    def deriv(self,dervar):
        raise NotImplementedError("The derivative of the poisson distribution is not implemented yet")
class Heaviside(Func):
    @property 
    def name(self):
        return "H"
    
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex,complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def hvs(num):
            try:
                return 1 if num>=0 else 0
            except TypeError:
                num=Complex(num)
                if num.isreal():
                    return hvs(num.real)
                else:
                    raise ValueError("Can' take sign of a complex number")
        return hvs
    @property
    def numpyfunc(self):
        def hvs(num):
            return np.heaviside(num,1)    
        return hvs
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        else:
            return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0*arg.deriv(dervar) #1/(arg**2+1)*arg.deriv(dervar) Ve skutecnosti by tu mela byt delta funkce  - implementovat?

class Max(Func):
    @property 
    def name(self):
        return "max"
    @property
    def numarg(self):
        return 2
    @property
    def evaluables(self):
        return int,float,Complex,complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def fmax(num1,num2):
            try:
                return max(num1,num2)
            except TypeError:
                num1,num2=Complex(num1),Complex(num2)
                if num1.isreal() and num2.isreal():
                    return fmax(num1.real,num2.real)
                else:
                    raise ValueError("Can't compare two complex numbers in function 'max'")
        return fmax
    @property
    def numpyfunc(self):
        def npmax(ar1,ar2):
            return np.maximum(ar1,ar2)
        return npmax
    def simplify(self):
        return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0*arg.deriv(dervar) #1/(arg**2+1)*arg.deriv(dervar) Ve skutecnosti by tu mela byt delta funkce  - implementovat?

class Min(Func):
    @property 
    def name(self):
        return "min"
    @property
    def numarg(self):
        return 2
    @property
    def evaluables(self):
        return int,float,Complex,complex
    @property 
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def fmin(num1,num2):
            try:
                return min(num1,num2)
            except TypeError:
                num1,num2=Complex(num1),Complex(num2)
                if num1.isreal() and num2.isreal():
                    return fmin(num1.real,num2.real)
                else:
                    raise ValueError("Can't compare two complex numbers in function 'min'")
        return fmin
    @property
    def numpyfunc(self):
        def npmin(ar1,ar2):
            return np.minimum(ar1,ar2)
        return npmin

    def simplify(self):
        return self
    def deriv(self,dervar):
        arg=self.attached[0]
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return 0*arg.deriv(dervar) #1/(arg**2+1)*arg.deriv(dervar) Ve skutecnosti by tu mela byt delta funkce  - implementovat?
def product(args): 
    
    return ft.reduce(mul,args,1)

def truncate_number(num,sig_nums=3):
    """Prevede cela cisla sama na sebe a desetinna cisla na cisla s desetinnou casti na sig_nums vyznamnych cislic"""
    inum=int(num)
    if inum==num: 
        return str(inum)
    else:
        base,decpart="{:.10f}".format(num).split(".")
        sig=decpart.strip("0")
        ind=decpart.find(sig)
        rounded=round(int(decpart[:ind+sig_nums+1]),-1)
        if rounded==10**(sig_nums+1):
            res=round(num,ind)
            ires=int(res)
            if ires==res: 
                return str(ires)
            else:
                return str(res)
        decimal=str(rounded).rstrip("0")
        #print(base,decpart,sig,int(decpart[:ind+sig_nums+1]),decimal,round(int(decpart[:ind+sig_nums+1]),-1))    
        
        return (base+"."+"0"*ind+decimal).rstrip(".")

def to_sigfigs(num,sig_nums=3): #  Vyznamne cislice s hezci vedeckou notaci
    import re
    s="{:.{sign}g}".format(num,sign=sig_nums)
    s=re.sub("e[\+\-](\d+)",lambda k: f"^{k[1]}",s)
    return s

def make_whole(x,y):
    #print("make whole:",x,y)
    orx,ory=x,y
    if x==0:
        return 0,y
    rx,ry=1/x,1/y
    #print("rx,ry",rx,ry)
    if int(rx)==rx:
        #print("int(rx)")
        p=rx*y
        if int(p)==p:
            #print("int(p)",p)
            return 1,p
        if int(ry)==ry:
            #print("int(ry)")
            return ry,rx
    if int(ry)==ry:
        p=x*ry
        if int(p)==p:
            return p,1

    #if abs(x)<abs(y):
    #    x,y=1/y,1/x
    for _ in range(5):
        if int(x)==x and int(y)==y:
            x,y=int(x),int(y)
            #print("looping",x,y)
            break
        x,y=x*10,y*10
    else:
        

        return orx,ory    
    #print("Returning",x,y)
    return x,y
    
class Monom(GeneralObject):
    @property
    def name(self):
        return "Monom"
    @property
    def evaluables(self):
        return int,float
    @property 
    def is_evaluable(self):
        return self.isconst
    def __init__(self,coef: float=0,vars :dict={}):

        self.isconst=True
        #print("Building Monom",coef,vars)
        if not vars:
            if isinstance(coef,Monom):
                coef,vars=coef.coef,coef.vars 
            elif isinstance(coef,GeneralObject):
                
                vars={repr(coef):(coef,1)}
                coef=1
        
        
        if coef==0:
            self.coef=0
            self.vars={}
            self.order=0
            self.varset=set()

        else:
            try:
                icoef=int(coef)
                if icoef==coef:
                    coef=icoef
            except TypeError:
                pass
            self.vars={}
            
            for name,obj_info in vars.items():
                if isinstance(obj_info,tuple):
                    #print("obj_info:",obj_info)
                    obj,order=obj_info
                else:
                    #obj,order=None,obj_info
                    obj,order=Variable(name),obj_info
                if order!=0:
                    self.vars[name]=(obj,order)
                    self.isconst=False
            self.coef=coef
            self.vars=dict(sorted(self.vars.items(),key=lambda item:ord(str(item[1][0])[0])-(100000 if isinstance(item[1][0],(Variable,Monom,Polynomial)) else 0)))
            self.order=sum(order for _,order in self.vars.values())
            self.varset={*self.vars.keys()}
        self.attached=[obj for obj,_ in self.vars.values()] 
         #   #print("From Monom:",coef,self.varset,self.order,self.isconst)
    
    def eval(self,*args,**kwargs)->float or int: #Povinna funkce, diky ktere Monom bude hrdou podtridou abstraktni tridy "Obecny vyraz" 
        res=self.coef
        
        for (obj,power) in self.vars.values():
        
            res*=obj.eval(**kwargs)**power#**power

        #print("returning Monom with coef,new_vars:",coef,new_vars)
        return res
     
        

    def at(self,**vals):
        new_coef=self.coef
        new_vars={**self.vars}
        for var,value in vals.items():
            if var in self.vars:
                new_coef*=value**self.vars[var][1]
                del new_vars[var]
        return Monom(new_coef,new_vars)
    @property
    def numpyfunc(self):
        def meval(*args,**kwargs):
            res=self.coef*1.0
            for arg,(_,power) in zip(args,self.vars.values()):
                res=res*arg**power#**power
            return res
        return meval   
    def issame(self,other):
        return self.vars==other.vars
    
    def inv(self):
        try:
            return Monom(1/self.coef if 1%self.coef else 1//self.coef,{
                var:-coef for var,coef in self.vars.items()})
        except ZeroDivisionError:
            #print("Deleni nulou!")
            raise
    def simplify(self):
        
        new_vars={}
        res=1
        for key,(obj,order) in self.vars.items():
            if isinstance(obj,Sqrt) and order%2==0:
                res*=obj.attached[0]**(order//2) #chybi simplify
            elif isinstance(obj,ObjPow) and order==1/obj.attached[1]:
                res*=obj.attached[0] #chybi simplify
            elif isinstance(obj,Func):
                res*=obj.simplify()**order
            else:
                    new_vars[key]=(obj,order)
        return Monom(self.coef,{var:(obj.simplify(),order) for var,(obj,order) in new_vars.items()})*res
    @ignore("Equation")
    def __add__(self,other):
        
        if not isinstance(other,Monom):
            if isinstance(other,(Polynomial,Quotient)):
                return other+self
            if not isinstance(other,GeneralObject):
                other=Monom(other,{})

            elif isinstance(other,Complex):
                return Complex(self+other.real,other.imag)
            elif isinstance(other,Vector):
                return Vector(*(self+comp for comp in other.values))
            elif isinstance(other,GeneralObject):
                other=Monom(other)

        if self.vars==other.vars:#self.issame(other):
            return Polynomial([Monom(self.coef+other.coef,self.vars)])
        else:
            return Polynomial([self,other])
            #raise ValueError("Can't add different monomials in Monom!")
    
    __radd__=__add__
    @ignore("Equation")
    def __sub__(self,other):
        return self+(-other)
    def __rsub__(self,other):
        return other+(-self)

    def monaddition(self,other):
        if self.vars==other.vars:#self.issame(other):
            return Monom(self.coef+other.coef,self.vars)
        else:
            raise ValueError("Monomy musi byt stejneho typu!")
    @ignore("Equation")
    def __mul__(self,other):
        #print("Multiplying",self,other)
        #print(type(other))
        if not isinstance(other,Monom):
            if isinstance(other,Variable):
                #updict={var:(obj,order+1 if other==obj else order) for var,(obj,order) in self.vars.items() }                
                #return Monom(self.coef,{var:(obj,order+1 if other==obj else order) for var,(obj,order) in self.vars.items() })
                return self*Monom(other)
            if isinstance(other,(Polynomial,Quotient)):
                return other*self
            if isinstance(other,(int,float)):
                coef=self.coef*other
                return Monom(coef,self.vars)
            if isinstance(other,Vector):
                return other*self
            if isinstance(other,Complex):
                return Complex(self*other.real,self*other.imag)
            if isinstance(other,(GeneralObject)):
                other=Monom(other)
        
        nv={**self.vars,**other.vars}
        nc=self.coef*other.coef
        if nc==0: 
            return Monom(0,{})
        else:
            for var1,(obj1,order1) in self.vars.items():
                for var2,(obj2,order2) in other.vars.items():
                    if var1==var2:
                        if obj1==obj2:
                            nv[var1]=(obj1,order1+order2)
                    #if var1=="const":
                    #   nv[var2]=order2
        res=Monom(nc,nv)
        #print("returning",res)
        return res
    def __hash__(self):
        return hash(repr(self))
    __rmul__=__mul__
    @ignore("Equation")
    def __truediv__(self,other):
        if isinstance(other, (int,float)):
           return Monom(self.coef/other if self.coef%other else self.coef//other,self.vars)
        return Polynomial(self)/other
    def __rtruediv__(self,other):
        return Quotient(other,self)
    def __neg__(self):
        return Monom(-self.coef,self.vars)

    def __pow__(self,exponent):
        if isinstance(exponent,Vector):
            return exponent.__class__(*(self**comp for comp in exponent))
        if type(exponent)!=int:
            return ObjPow(self,exponent)
            #raise TypeError("Exponent musi byt cele cislo!")
        if exponent>0:
            return Monom(self.coef**exponent,{var:(obj,order*exponent) for var,(obj,order) in self.vars.items()})            
        elif exponent<0:
            return Quotient(1,self**(-exponent))
        else:
            return Monom(1)

    def deriv(self,dervar):
        res=0 
        for var,(obj,order) in self.vars.items():
            rest=Monom(self.coef,{nvar:res for nvar,res in self.vars.items() if nvar!=var})
            
            dinner=obj.deriv(dervar)
            if order-1!=0:
                douter=Monom(order,{var:(obj,order-1)})
            else:
                douter=order
            res+=douter*dinner*rest
        
        return res
    @property
    def partials(self):
        
        res=[]
        for var,(obj,order) in self.vars.items():
            rest=Monom(self.coef,{nvar:res for nvar,res in self.vars.items() if nvar!=var})
           # print("rest",rest)
            if order-1!=0:
                douter=Monom(order,{var:(obj,order-1)})
               # print("douter",douter)
            else:
                douter=order
            res.append(douter*rest)
        
        return res

    def divgrad(self):
        return Polynomial(self).divgrad()

    def __eq__(self,other):

        if isinstance(other,(float,int)):
            if self.isconst: 
                return self.coef==other
            else:
                return False
        if isinstance(other,Complex) and other.isreal():
                return self==other.real
        if isinstance(other,Monom):
            return self.coef==other.coef and self.vars==other.vars
        if isinstance(other,(Polynomial,Quotient)):
            return other.__eq__(self)
        return super().__eq__(other)
    
    def __repr__(self):
        #return f"{module_prefix}Monom({self.coef!r},{self.vars!r})"
        #return f"{module_prefix}M({self.coef!r},{self.vars!r})"
        return f"{module_prefix}M({self.coef!r},{dict((key,order) for key,(_,order) in self.vars.items())!r})"
        #self.vars=(("x",1),("y",5))
    def __int__(self):
        if self.isconst:
            return int(self.coef)
        else:
            raise ValueError("Can't convert a non-constant monomial to int")
    def __float__(self):
        if self.isconst:
            return float(self.coef)
        else:
            raise ValueError("Can't convert a non-constant monomial to float")
    @classmethod
    def fromrepr(cls,text):
        if not text.startswith("M("):
            raise ValueError("Wrong format! The string must start with \"M(\" ")
        else:
            
            return eval("Monom("+text[2:])
        
    def __str__(self):
        s=""
        if self.isconst:
            return truncate_number(self.coef)
        else:
            if self.coef>=0:
                pref=""
            else:
                pref="-"
            #print("#printing self...",*(str(key)+":"+str(obj)+"..."+str(order) for key,(obj,order) in self.vars.items()))
            trunc_coef=truncate_number(abs(self.coef))
            s+=f"{pref}{trunc_coef if trunc_coef!='1' else ''}"
            for _,(obj,order) in self.vars.items():
                    s+=str(obj)
                    if order!=1:
                        if order>0:
                            s+=f"^{order}"
                        else:
                            s+=f"^({order})"
        
            return(s)
        
           

class Polynomial(GeneralObject):
    @property
    def name(self):
        return "Polynomial"
    @property
    def evaluables(self):
        return int,float
    @property
    def is_evaluable(self):
        return all(mon.is_evaluable for mon in self.monlist)

    def __init__(self,monlist):
        
        nl=[]
        if isinstance(monlist,Polynomial):
            monlist=monlist.monlist
        try:
            iter(monlist)
            if isinstance(monlist,Complex):
                
                if monlist.isreal():
                    monlist=[monlist.real]
                elif monlist.isimag():
                    monlist=[Complex(0,1)*monlist.imag]
                else:
                    monlist=[monlist]
                #print("monlist:",monlist)
            elif isinstance(monlist,Vector):
                monlist=[monlist]
        except TypeError:
            if isinstance(monlist,Polynomial):
                monlist=monlist.monlist
            elif isinstance(monlist,Quotient):
                if monlist.q.isconst:
                    const=(monlist.q.monlist)[0].coef
                    monlist=[mon/const for mon in monlist.p.monlist]

                else:
                    raise ValueError("Can't convert quotient with non-constant denominator to a polynomial!")
            else:
                
                monlist=[monlist]
                                       
        monlist=list(monlist)
        
        for mon in monlist:
            if not isinstance(mon,Monom):
                mon=Monom(mon)
            same=False
            
            for ind,savedmon in enumerate(list(nl)):
                if mon.vars==savedmon.vars:#mon.issame(savedmon):
                    
                    same=True
                    res=mon.monaddition(nl[ind])
                    
                    if res==0:
                        del nl[ind]
                    else:
                        nl[ind]=res
            if not same:
                if mon!=0:
                    nl.append(mon)
        
        self.monlist=sorted(nl,key=self.msum)
        if not self.monlist:
            self.monlist.append(Monom(0,{}))
        
        self.order=max(mon.order for mon in self.monlist)
        self.varset={var for mon in self.monlist for var in mon.vars.keys()}
        self.isconst=all(mon.isconst for mon in self.monlist) 
        self.attached=self.monlist
        self.obj_dict={}
        for var in self.varset:
            #print("var:",var)
            for mon in self.monlist:
                #print("mon:",mon)
                if var in mon.vars:
                    self.obj_dict[var]=mon.vars[var][0]
                    break

    def eval(self,*args,**kwargs):
        if len(self.varset)!=1:
            return sum(mon.eval(*args,**kwargs) for mon in self.monlist)
        else:
            return self.horneval(*args,**kwargs)

    def horneval(self,*args,**kwargs): #Horner's scheme
        var=next(iter(self.varset))
        for mon in self.monlist:
            if not mon.isconst:
                val=mon.vars[var][0].eval(*args,**kwargs)
                break
        
        an=self.get_coefs(var)
        res=an[-1]
        for ai in reversed(an[:-1]):
            res=res*val+ai
        return res

    @property
    def numpyfunc(self):
        def meval(*args):
            res=0.0
        
            for arg,_ in zip(args,(mon for mon in self.monlist)):
                #res+=arg
                res=res+arg
            return res
        return meval
    def simplify(self): 
        sinlike=Sin,Sinh;coslike=Cos,Cosh;signs=(1,-1)
        expr=Polynomial(self)
        for S,C,sign in zip(sinlike,coslike,signs):
            pref,rest=expr.factor_variables()
            found=False # sin^2+cos^2=1
            for mon in rest.monlist:
                for (obj,power) in mon.vars.values():
                    if type(obj) is C and power>=2:
                        monmark=mon
                        newmon=mon/obj**2
                        found=True
                        break
                        #coef=monmark.coef#;d={name2:(obj2,power2) for name2,(obj2,power2) in  mon.vars.items() if obj2!=obj };
                        #arg=obj.attached[0]
                        
                    
            if found : 
                found=False
                for mon in rest.monlist:
                    for obj,power in mon.vars.values():
                        if type(obj) is S and power>=2:
                            #d2={obj2:power2 for obj2,power2 in  mon.vars.values() if obj2!=obj }
                            if newmon==sign*mon/obj**2:
                                found=True
                                monmark2=mon
                        #if d==d2:
                        #    monmark2=mon
            if found:
                new_monlist=[mon for mon in rest.monlist if mon not in [monmark,monmark2]]+[newmon]
                expr=Polynomial(pref.simplify()*sum(mon.simplify() for mon in new_monlist)       )
                
        return sum(mon.simplify() for mon in expr.monlist)        
                 
    def factor_variables(self):
        #print("Factoring")
        d=self.monlist[0].coef
        #print(self.varset)
        
        if self==0:
            return (Monom(0),Polynomial(0))
        for mon in self.monlist[1:]:
           # print("entering:",d,mon.coef)
            d=gcd(d,mon.coef,tol=self.tol)
            #print("leaving:",d)
        #print("d=",d,self.monlist,",",self.monlist[1:])
        pref_vars={var:(self.obj_dict[var],min(mon.vars.get(var,(0,0))[1] for mon in self.monlist)) for var in self.varset}
        prefactor=Monom(d,pref_vars)
        new_coefs=(mon.coef//d  if d!=1 else mon.coef for mon in self.monlist)
        new_vars=({var:(obj,order-pref_vars[var][1]) for var,(obj,order) in mon.vars.items() if order-pref_vars[var][1]!=0}  for mon in self.monlist)
        rest=Polynomial(Monom(coef,var) for coef,var in zip(new_coefs,new_vars))

        #print("Returning from factor","|".join(str(el) for el in (prefactor,rest,prefactor.coef,prefactor.isconst)))
        return (prefactor,rest)
   
    @classmethod
    def convert(cls,obj):
        if isinstance(obj,Polynomial):
            return obj
        if isinstance(obj,(float,int)):
            #print(obj)
            return SPolynomial("",obj)
        elif isinstance(obj,Monom):
                return Polynomial([obj])
        else:
            raise ValueError(f"Can't convert object {obj} with type {type(obj)} to Polynomial")
    def msum(self,mon):
    
    #print(arglist,len(arglist))
 
        return -(mon.order*100000-len(mon.vars)*1000-0.01*sum(ord(str(key[0])) for key in mon.vars.keys()))

    def at(self,**vars):
        return sum(mon.at(**vars) for mon in self.monlist)

    @ignore("Equation")
    def __add__(self,other):
        #print("adding",self,other,type(self),type(other),self.monlist,other.monlist)
        
        if isinstance(other,Polynomial):
            return Polynomial(self.monlist+other.monlist)
        elif isinstance(other,Monom):
            return Polynomial(self.monlist+[other])    
        elif isinstance(other,Quotient):    
            return other+self
        elif isinstance(other,Complex):
            return Complex(self+other.real,other.imag)
        elif isinstance(other,Vector):
            return Vector(*(self+comp for comp in other.values))
        else:
            return self+Monom(other)
    
    __radd__=__add__

    @ignore("Equation")
    def __mul__(self,other):
        
        if isinstance(other,Polynomial):
            return Polynomial(mon1*mon2 for mon1 in self.monlist for mon2 in other.monlist)

        if isinstance(other,(int,float)): #nebo cokoliv jineho 
            return Polynomial(other*mon for mon in self.monlist)
        if isinstance(other,Monom):
            return Polynomial(mon1*mon2 for mon1 in self.monlist for mon2 in [other])    
           # return Polynomial(mon1*other for mon1 in self.monlist)    
        if isinstance(other,Quotient):
            return other*self
        if isinstance(other,(Complex,Vector)):
            return other.__class__(*(self*comp for comp in other.values))
        return self*Monom(other)
         
                
    __rmul__=__mul__

    def __neg__(self):
        nl=[(-1)*mon for mon in self.monlist]
        return Polynomial(nl)

    @ignore("Equation")	
    def __pow__(self,exponent):
        if isinstance(exponent,Vector):
            return exponent.__class__(*(self**comp for comp in exponent))
        if not isinstance(exponent,int): #or exponent<0:
            return ObjPow(self,exponent)
            #raise TypeError("Exponent musi byt cele cislo!")
        if exponent==0:
            return Polynomial(1)
        if len(self.monlist)==1:
            if exponent>0:
                return Polynomial(self.monlist[0]**exponent)
            else:
                return Quotient(1,self**(-exponent))
            
        else:
            if exponent>0:
                return ft.reduce(mul,(self for _ in range(exponent)))
            else:
                return Quotient(1,self**(-exponent))

    @ignore("Equation")       
    def __sub__(self,other):
        return self+(-other)
    def __rsub__(self,other):
        return other+(-self)

    @ignore("Equation")
    def __truediv__(self,other):
        #print("Pol div",repr(self),"|",repr(other))
        if isinstance(other,(int,float)):
            if other==0:
                raise ZeroDivisionError("Deleni nulou!")
            else:
                return self*(1/other if 1%other else 1//other)
        if other==self:
            
            
            return Quotient(1)
        
        if isinstance(other,Monom):
            #return self*other.inv()
            return self/Polynomial(other)
        if isinstance(other,Quotient):
            #print("other.inv",other.inv())
            return other.inv()*self
        #if isinstance(other,GeneralObject):
        #    other=Monom(other)    
        #print(type(self),type(other))        
        #print("just passing through...")
        return Quotient(self,other) 

    def __rtruediv__(self,other):

        return Quotient(other)/self
    
    def derivo(self,var):
        nl=[]
        for mon in self.monlist:
            #print(var,mon.vars,mon.coef)
            if var not in mon.vars:
                continue
            new_vars={**mon.vars}
            new_vars[var]-=1
            #print(mon.coef,mon.vars[var]) 
            nl.append(Monom(mon.coef*mon.vars[var][1],new_vars))
        return Polynomial(nl)

    def deriv(self,dervar):
        nl=Polynomial(0)
        for mon in self.monlist:
            nl+=mon.deriv(dervar)

        return nl
    def divgrad(self):
        nl=[]
        for mon in self.monlist:
            if mon.isconst:
                continue
            for var,(_,order) in mon.vars.items():
                if order==0:
                    continue
                new_vars={**mon.vars}
                new_vars[var]-=1
                nl.append(Monom(mon.coef*order,new_vars))
        return Polynomial(nl)
    
    def __repr__(self):
        #return(f"{module_prefix}Polynomial(["+",".join((f"{mon!r}" for mon in self.monlist))+"])")
        return(f"{module_prefix}P(["+",".join((f"{mon!r}" for mon in self.monlist))+"])")

    def __eq__(self,other):
        #print(self.order)
        if isinstance(other,Complex):
            if other.isreal():
                return self==other.real
        if self.isconst:
            if isinstance(other,(float,int)):
                #print(self.monlist[0],other)
                return self.monlist[0]==other
        
        if isinstance(other,Monom):
            return self==Polynomial(other)
        
        if isinstance(other,Polynomial):
            #print(set(self.monlist))
            #print(set(other.monlist))#,"comp")
            #print(self.monlist[0].coef==other.monlist[0].coef)
            #print(self.monlist[0].vars==other.monlist[0].vars)
            #for ((key,(obj,order)),(key2,(obj2,order2))) in zip(self.monlist[0].vars.items(),other.monlist[0].vars.items()):
             #   #print(key,key2,key==key2)
             #   #print(obj,obj2,repr(obj),repr(obj2),obj==obj2)
              #  #print(order,order2,order==order2)
            #print(self.monlist[0].vars==other.monlist[0].vars)
            return set(self.monlist)==set(other.monlist)
        elif isinstance(other,Quotient):
            return other==self
        return super().__eq__(other)

    def __str__(self):
        #smonlist=sorted(self.monlist,key=lambda mon: product((order for order in mon.vars.values())))
        s=""
        for mon in self.monlist:
            if mon.coef!=0:
                if s=="":
                    s+=str(mon)
                    
                else:
                    if mon.coef>0:
                        s+=" + "+str(mon)
                    else:
                        s+=" - "+str(mon)[1:]   
        if s=="":
            s="0"               
        return s

    
    def __int__(self):
        for mon in self.monlist:
            if  mon.isconst:
                return int(mon.coef)
        raise ValueError("Can't convert a non-constant polynomial to int")

    def __float__(self):
        for mon in self.monlist:
            if  mon.isconst:
                return float(mon.coef)
        raise ValueError("Can't convert a non-constant polynomial to float")

    @classmethod
    def fromrepr(cls,text):
        return eval(text.replace("P","Polynomial").replace("M","Monom"))
    def get_coefs(self,var,at_least=None):
       # #print(self.varset,self.order)

        max_order=max(at_least,self.order) if at_least!=None else self.order
        #coefs={order:0 for order in range(max_order+1)}
        coefs=[0 for _ in range(max_order+1)]
        for mon in self.monlist:
            if mon.vars=={}:
                coefs[0]=mon.coef
            else:
                coefs[mon.vars[var][1]]=mon.coef
        return coefs

    def get_linear_coefs(self,varset=None,at_least=None):
        
        
        constcoef=0
        if varset is None:
            varset=self.varset
       
        coefdict={var:0 for var in varset}
        for mon in self.monlist:
            if mon.vars=={}:
                constcoef=mon.coef
            else:
                for var in varset:
                    if var in mon.vars:
                        coefdict[var]=mon.coef
                        break
        return coefdict,constcoef
            
    def leading_term(self,var): 
            #return (0,pol.monlist[0].coef) if pol.isconst else max(filter(lambda d: d[1]!=0,(d for d in pol.get_coefs(var).items()) ))                    
            #return (0,pol.monlist[0].coef) if pol.isconst else max(filter(lambda d: d!=0,(d for d in pol.get_coefs(var))))
            #print(self.get_coefs(var))
            #input()
            for order,coef in enumerate(self.get_coefs(var)[::-1]):
                if coef!=0:
                    #print("returning",self.order-order,coef)
                    return self.order-order,coef
            return 0,0

    def longdiv(self,other):
        #print("longdiving with",(self),(other))
        if isinstance(other, (int,float,Monom,Variable)):
            other=Polynomial(other)
        elif isinstance(other,Quotient):
            if other.q.isconst:
                other=Polynomial(other.p/other.q.monlist[0])
            else:
                raise TypeError("Can't use the long division algorithm to divide by a quotient")
        if not self.varset:
            if not other.varset:
                num=self.monlist[0].coef
                den=other.monlist[0].coef
                return Polynomial(num//den),Polynomial(num%den),[(self,0),(num%den,num//den)]
            else:
                return Polynomial(0),self,[(self,0),(self,0)]
        if len(self.varset)>1:
                raise NotImplementedError(r"The polynomials must be of the same type and in one variable,"
        r"the more general case hasn't been implemented yet")
        if other.varset:
            if  self.varset!=other.varset:
                raise NotImplementedError(r"The polynomials must be of the same type and in one variable,"
        r"the more general case hasn't been implemented yet")
        
        var=next(iter(self.varset))

        p=self
        q=other            

        s=SPolynomial(var,0)
        partial_results=[(p,0)]

        porder,pcoef=p.leading_term(var)
        #print(porder,pcoef)
        qorder,qcoef=q.leading_term(var)
        
        while porder>=qorder:
            #print(p,":",q) 
            #print(porder,pcoef,qorder,qcoef)
            #input()
            #print(p,pcoef,q,qcoef,q.obj_dict)    
            if p.obj_dict:
                newterm=Monom(pcoef/qcoef if pcoef%qcoef else pcoef//qcoef,{var:(p.obj_dict[var],porder-qorder)})
            else:
                newterm=Monom(pcoef/qcoef if pcoef%qcoef else pcoef//qcoef)
            #print(pcoef,qcoef,"...",s,p,newterm*q,type(newterm),type(q),type(newterm*q))
            #print(" | ".join(str(el) for el in (var,repr(p),q,repr(newterm),q*newterm,newterm*q)))
            #input()
            #print(newterm,p,s)
            s+=newterm
            p-=newterm*q
            #print("->",p,s)
            #print(newterm,newterm*q)
            #print("s,r:",s,"|",p)
            #print(p,p.order)
            #print(repr(p),p==0,str(p))
            partial_results.append((p,s))
            if p==0:
                break
            
            porder,pcoef=p.leading_term(var)
            
         #   if self.order<other.order:
          #      return SPolynomial("",0.0),self
        #print(self.varset)
        
        if len(partial_results)==1:
            partial_results.append((self,0))
        return s,p,partial_results
       # except AttributeError as e:
        #    #print(e)
         #   #print("The argument must a monomial, a polynomial or a polynomial quotient!")

class SMonom(Monom): 
    def __init__(self,coef,symb,order):
        super().__init__(coef,{symb:order})

class SPolynomial(Polynomial):
    def __init__(self,symb,*coefs):
        var=Variable(symb)
        super().__init__((Monom(coef,{symb:(var,i)}) for i,coef in enumerate(coefs)))

class Quotient(GeneralObject):
    @property
    def name(self):
        return "Quotient"
    @property
    def evaluables(self):
        return int,float
    @property
    def is_evaluable(self):
        return self.p.is_evaluable and self.q.is_evaluable

    def __init__(self,p:Polynomial,q:Polynomial=None,cancel_common=True):
        
        if p is None:
            raise ValueError("P can't be None!")
        if q is None:
            q=1
        #if isinstance(p,Quotient) and isinstance(q,Quotient):
        #    p,q=p.p*q.q,p.q*q.p
        #elif isinstance(p,Quotient):
        #    p,q=p.p,p.q*q
        #elif isinstance(q,Quotient):
        #    p,q=p*q.q,q.p
        
        if isinstance(p,Quotient):
            if isinstance(q,Quotient):
                p,q=p.p*q.q,p.q*q.p
            else:
                p,q=p.p,p.q*q
        elif isinstance(q,Quotient):
            p,q=p*q.q,q.p
        if q==0:
            raise ZeroDivisionError("Division by zero!")   
        
        if p==0:
            p=Polynomial(0)
            q=Polynomial(1)
        else:
                
            if q is None:
                q=SPolynomial("",1)
            
            #print(p,q)
            if not isinstance(p,Polynomial): 
                p=Polynomial(p)
            
            if not isinstance(q,Polynomial):
                q=Polynomial(q)
            
            if p==q:
                
                p=Polynomial(1)
                q=Polynomial(1)
            elif p.isconst and q.isconst:
                #print("are const")
                a=p.monlist[0].coef
                b=q.monlist[0].coef
                #print("a,b",a,b)
                a,b=make_whole(a,b)
                #print("a,b",a,b)
                com_div=gcd(a,b)
                a,b=a/com_div,b/com_div
                p=Polynomial(a)
                q=Polynomial(b)
          
            else:
                #print("We are here",p,",",q)
                #print("middling quotient with...",repr(p)," | ",repr(q)," | ")                      
                if p.varset==q.varset and len(p.varset)==1:
                 #    #print("longdiving",p,q)
                        

                    if p.order>=q.order:
                        s,r,_=p.longdiv(q)                
                        if r==0:
                            p=s
                            q=Polynomial(1)
                    else:
                        s,r,_=q.longdiv(p)                
                        if r==0:
                            p=Polynomial(1)
                            q=s
                    #print(s,r)
                 #   p,q=self.cancel_and_divide(p,q)

                else:
                    
                    if  cancel_common:
                        #print("cancel_common")    
                        p,q=self.cancel_and_divide(p,q)
                        #print("After cancelling, got",p,q)
        if q.monlist[0].coef<0:
            self.p,self.q=-p,-q            
        else:
            self.p,self.q=p,q
            
        self.varset=self.p.varset.union(self.q.varset)
        self.isconst=self.p.isconst and self.q.isconst
        self.attached=[self.p,self.q]
        #print("exiting quotient with...",repr(self.p)," | ",repr(self.q))                      
        #print("Finished Quotient with ",self.p,self.q)
            #print("Built Quotient with p,q:,type(p),type(q),varset,isconst:",self.p,self.q,type(p),type(q),self.varset,self.isconst)

    def cancel_and_divide(self,p,q):
        #print("cancelling",p,"|",q)
        ppref,prest=p.factor_variables()
        #print("factoring q",p,"|",q)
        qpref,qrest=q.factor_variables()
       # #print("|".join(str(el) for el in (ppref,qpref,prest,qrest)))
        #input()
        
        #print(" | ".join(map(str,(ppref,qpref,prest,qrest))))
        pcoef,qcoef=ppref.coef,qpref.coef
        #print(f"d=gcd({pcoef},{qcoef})")
        d=gcd(pcoef,qcoef)

        pcoef,qcoef=pcoef/d,qcoef/d #pozor, co to udělá když jedno z toho bude float?
        obj_dict={**p.obj_dict,**q.obj_dict}
        common_vars={var:(obj_dict[var],(min(mon.vars.get(var,(None,0))[1] for mon in (ppref,qpref)))) for var in ppref.varset.union(qpref.varset)}
        newpvars,newqvars=({var:(obj,order-common_vars[var][1]) for var,(obj,order) in mon.vars.items() if order-common_vars[var][1]!=0} for mon in (ppref,qpref))
        ppref,qpref=Monom(pcoef,newpvars),Monom(qcoef,newqvars)
        
        #print("v1:",prest,qrest)
        #input()
       # #print(all(mon1.coef==mon2.coef for mon1,mon2 in zip(prest.monlist,qrest.monlist)))
      #  #print(prest==qrest,repr(prest),repr(qrest))
        
       # #print(v1)
        #v2=v1*ppref
        #print("v2 pass",v2,type(v2))
        #v3=v2/qpref
        #print("v3 pass",v3,type(v3))
        #print("Sending",prest,qrest,"and ",ppref,qpref," to newz",pcoef,qcoef)

        z1=Quotient(ppref,qpref,cancel_common=False)
        z2=Quotient(prest,qrest,cancel_common=False)
        
        newz=Quotient(z1.p*z2.p,z1.q*z2.q,cancel_common=False)
        return newz.p,newz.q

    def eval(self,*args,**kwargs):
        #print("---",repr(self.p.eval(*args,**kwargs)),"/",repr(self.q.eval(*args,**kwargs)))
        #print(repr(self.p.eval(*args,**kwargs)))
        #print(repr(self.q.eval(*args,**kwargs)))
        return self.p.eval(*args,**kwargs)/self.q.eval(*args,**kwargs)
    @property
    def numpyfunc(self):
        return truediv
    def at(self,**vars):
        res=self.p.at(**vars)/self.q.at(**vars)
        if isinstance(res,float):
            if int(res)==res:
                return int(res)
        else:
            return res
    @ignore("Equation")
    def __add__(self,other):
     #   #print("type other",type(other),isinstance(other,Quotient))
        if not isinstance(other,Quotient):
            if isinstance(other,Complex):
                return Complex(self+other.real,other.imag)
            if isinstance(other,Vector):
                return other.__class__(*(self+comp for comp in other.values))
          #  if isinstance(other,Monom):
          #      other=Polynomial(other)
            else: #isinstance(other,(Polynomial,int,float)):
                return Quotient(self.p+other*self.q,self.q)    
        if self.q==other.q:
            return Quotient(self.p+other.p,self.q)    
        else:
            return Quotient(self.p*other.q+self.q*other.p,self.q*other.q)
    __radd__=__add__
   

    @ignore("Equation")
    def __sub__(self,other):
        return self+(-other)

    def __rsub__(self,other):
        return other+(-self)
    @ignore("Equation")
    def __mul__(self,other):
        #if isinstance(other,(Variable,Polynomial,Monom,int,float)):# or type(other) in [int,float]:
        if not isinstance(other,Quotient):# or type(other) in [int,float]:
            #return Quotient(self.p*other,self.q)
            if isinstance(other,(Complex,Vector)):
                return other.__class__(*(self*comp for comp in other.values))
            else:
                other=Quotient(other)
        if self.p!=0:
            
            if self.p==other.q and self.q==other.p:
                #print("1545")   
                return Quotient(1)
            elif self.p==other.q:
                return Quotient(other.p,self.q)
            elif other.p==self.q:
                #print("1550")   
                return Quotient(self.p,other.q)
            elif self.p*other.p==self.q*other.q:
                return Quotient(1)
            else: 
                #print("here")   
                return Quotient(self.p*other.p,self.q*other.q)
        else:
            return Quotient(0)
    __rmul__=__mul__

    def __pow__(self,exponent):
        if isinstance(exponent,Vector):
            return exponent.__class__(*(self**comp for comp in exponent))
        if not isinstance(exponent,(int,float)):
            return ObjPow(self,exponent)
        if exponent>0:
            return Quotient(self.p**exponent,self.q**exponent)
        elif exponent<0:
            if self.p==0:
                raise ZeroDivisionError("Nula nema zaporne mocniny!")
            else:
                return Quotient(self.q**(-exponent),self.p**(-exponent))
        else: #exponent==0
            return Quotient(0)

    

    def __repr__(self):

        cs=self
        #return f"{module_prefix}Quotient({cs.p!r},{cs.q!r})"
        return f"{module_prefix}Q({cs.p!r},{cs.q!r})"

    def __str__(self):
        cs=self
        if self.q.isconst and self.q.monlist[0].coef<0:
                p,q=-self.p,-self.q
        else:
            p,q=self.p,self.q
        
        p,q=cs.p,cs.q
        pstr=str(p)
        qstr=str(q)
        
        pstrlen=max(len(pstr.split(" + ")),len(pstr.split(" - ")))
        upar="()" if pstrlen>1 or len(self.p.varset)>1 else ["",""]
        if q.isconst or (len(q.monlist)==1  and len(q.monlist[0].vars)==1) and q.monlist[0].coef==1: 
                lpar=["",""]
        else:
            lpar="()"
        #qstrlen=max(len(qstr.split(" + ")),len(qstr.split(" - ")))
        #print("debug",self.p,self.q)
        if q==1:
            return pstr
        
        return ("{upar[0]}{up}{upar[1]}/{lpar[0]}{down}{lpar[1]}".format(
            up=pstr,down=qstr,upar=upar,lpar=lpar))
    def __format__(self,fmt):
        if fmt=="q":
            pstr=str(self.p)
            qstr=str(self.q)

            pstrlen=max(len(pstr.split(" + ")),len(pstr.split(" - ")))
            qstrlen=max(len(qstr.split(" + ")),len(qstr.split(" - ")))
            upar="()" if pstrlen>1 else ["",""]
            lpar="()" if qstrlen>1 else ["",""]
        #print("debug",self.p,self.q)
            return ("{upar[0]}{up}{upar[1]}/{lpar[0]}{down}{lpar[1]}".format(
            up=pstr,down=qstr,upar=upar,lpar=lpar))
        else: 
            return str(self)
    def __int__(self):
        if self.isconst:
            return int(float(self.p)/float(self.q))
        else:
            raise ValueError("Can't convert a non-constant quotient to int")
    def __float__(self):
        
        if self.isconst:
            return float(float(self.p)/float(self.q))
        else:
            raise ValueError("Can't convert a non-constant quotient to float")
        
    def inv(self):
        #print(f"self {self},|{self.p},{self.q} | invself {Quotient(self.q,self.p)}")
        #print("from inv:",self.p,self.q)
        return Quotient(self.q,self.p,cancel_common=False)
    def simplify(self):
        return self.p.simplify()/self.q.simplify()
    def __neg__(self):
        return Quotient(-self.p,self.q)
    def deriv(self,var):
        return (self.p.deriv(var)*self.q-self.p*self.q.deriv(var))/(self.q**2)
    @ignore("Equation")
    def __truediv__(self,other):
        #print("Quot div",self,"|",other)
        

        if not isinstance(other,Quotient):
            #print("self,other",self,other)
            other=Polynomial(other)

        if other==0:
            raise ZeroDivisionError("Deleni nulou!")
        else: 
            #print("here",self,other,type(other))#,type(other.inv()),self*other.inv())
            return Quotient(self.p,self.q*other)#other.inv()

    def __rtruediv__(self,other):
        
        if self.p==0:
            raise ZeroDivisionError("Deleni nulou!")
        else:
            if not isinstance(other,Quotient):
                other=Polynomial(other)
            return other/self
    def __eq__(self,other):
        if  not isinstance(other,Quotient):
            #print(self,other)
            #print(",".join(map(repr,(self.p,self.q,other,other*self.q))))
            if isinstance(other,Complex) and other.isreal():
                return self==other.real
                
            if self.q==1:
                return self.p==other
            else:
                if isinstance(other,(int,float,complex,GeneralObject)):    
                    return self.p==other*self.q
                else:
                    return False
        else:
            return self.p*other.q==self.q*other.p
    @property                        
    def real(self):
        p,q=self.p,self.q
        pr=p.real if isinstance(p,Complex) else p;qr=q.real if isinstance(q,Complex) else q
        pi=p.imag if isinstance(p,Complex) else 0;qi=q.imag if isinstance(q,Complex) else 0
        
        return Quotient(pr*qr+pi*qi,qr**2+qi**2)
    @property
    def imag(self):
        p,q=self.p,self.q
        pr=p.real if isinstance(p,Complex) else p;qr=q.real if isinstance(q,Complex) else q
        pi=p.imag if isinstance(p,Complex) else 0;qi=q.imag if isinstance(q,Complex) else 0
        return Quotient(pi*qr-pr*qi,qr**2+qi**2)

class Vector(GeneralObject):
    @property
    def name(self):
        return "Vector"
    @property
    def evaluables(self):
        return (int,float)
    def is_evaluable(self):
        return True
    def eval(self,*args,**kwargs):
        return self.__class__( *(comp.eval(*args,**kwargs) if isinstance(comp,GeneralObject) else comp for comp in self.values))
    @property
    def numpyfunc(self):
        def arize(*args):
            #print(*(np.shape(arg) for arg in args))
            return np.array(args)
        return arize
        
    def __init__(self,*values):
        #print("Vector got values",values)
        self.attached=self.values=list(values)
        self.dim=len(self.values)
        #self.tol=1.e-10
    
    def __len__(self):
        return len(self.values)
    @ignore("Equation")#,"RandomSample")
    def __add__(self,other):
        if isinstance(other,Vector): #uncomplex here
            if self.dim!=other.dim:
                raise ValueError("The vectors must have same dimension!")
            return Vector(*(val1+val2 for val1,val2 in zip(self.values,other.values)))
        else:
            return Vector(*(val+other for val in self.values) )
    __radd__=__add__

    def __neg__(self):
        return Vector(*(-val for val in self.values))
    @ignore("Equation")
    def __sub__(self,other):
        return self+(-other)
    def __rsub__(self,other):
        return -self+other
    @ignore("Equation")
    def __mul__(self,other):
        #if isinstance(other,(int,float)):
         #   raise NotImplementedError("A vector can only be multiplied by a number! The Lie algebra is not implemented yet")
        #else:
        #if isinstance(other,(Complex,Vector)):
        #    return other.__class__(*(self*val for val in other) )
        if isinstance(other,Vector):
            return Vector(*(self*val for val in other) )
        elif isinstance(other,Complex):
            return Vector(*(Complex(comp*other.real,comp*other.imag) for comp in self) )
        else:
            return Vector (*(other*val for val in self) )

    __rmul__=__mul__
    @ignore("Equation")
    def __truediv__(self,other):
        if isinstance(other,Vector):
            raise ValueError("Can't divide two vectors!")
        else:
            return Vector(*(comp/other for comp in self.values))
    def __pow__(self,exponent):
        if isinstance(exponent,Vector):
            return exponent.__class__(*(self**comp for comp in exponent))
        """if isinstance(exponent,(int,float)) and  int(exponent)==exponent:
            if exponent>0:
                v=self
                for _ in range(exponent-1):
                    v=self*v
                return v
            elif exponent==0:
                return Monom(1)
        else:
            raise ValueError("The exponent must be an integer")
        """
        return self.__class__(*(comp**exponent for comp in self))
    def __rpow__(self,base):
        return self.__class__(*(base**comp for comp in self))

    def __rtruediv__(self,other):
        return Vector(*(other/comp for comp in self.values))
        
    def norm(self):
        return Sqrt(self.dot(self)).eval()
    
            

    def T(self):
        subvecs=[]
        try:
            l=self.values[0].dim
            check=all(comp.dim==l for comp in self.values)
            if not check:
               raise ValueError("Can't transpose object with unequal number of elements in each column!")
            veclen=self.values[0].dim
            for ind in range(veclen):
                subvecs.append(Vector(*(vec.values[ind] for vec in self.values)))
            return Vector(*subvecs)
        except (TypeError,AttributeError): 
            return self
        
    def cross(self,other):
        if self.dim==2:
            vec1=Vector(*self.values,0)
        else:
            vec1=self
        if other.dim==2:
            vec2=Vector(*other.values,0)
        else:
            vec2=other
        #print(vec1.dim,vec1,vec2.dim,vec2)
        if vec1.dim==vec2.dim!=3:
            raise ValueError("The vector product of two vectors satisfying Jacobi's identity is only defined in dimension 3!")
        x1,y1,z1=vec1.values
        x2,y2,z2=vec2.values
        return Vector(y1*z2-y2*z1,z1*x2-z2*x1,x1*y2-x2*y1)

    def dot(self,other):
        if not isinstance(other,Vector):
            raise ValueError("The right argument has to be a vector!")
        if self.dim!=other.dim:
            raise ValueError("The vectors need to have identical dimensions!")
        return sum(val1*val2 for val1,val2 in zip(self.values,other.values))
    def matmul(self,obj2):
            try:
                #print(obj1,obj2)
                return self.dot(obj2.T())
            except AttributeError:
                raise TypeError("The arguments must be vectors!")
    def __repr__(self):
        return f"Vector("+",".join(repr(value) for value in self.values)+")"
    def __str__(self):
        #print("vvals:",self.values)
        #if len(self.values)==1:
        #    return f"V({self.values[0]})"
        #else:
            return f"( "+" , ".join(truncate_number(value) if isinstance(value,float) else "".join(map(str,value))
        if isinstance(value,tuple) else str(value) for value in self.values)+" )"
    def __eq__(self,other):
        if not isinstance(other,Vector):
            return  all(comp==other for comp in self)
        #print("Comparing",repr(self),"and",repr(other))
        #print(self.values,other.values,self.values[0]-other.values[0],self.values[1]-other.values[1])
        for c1,c2 in zip(self.values,other.values):
            if c1!=c2:
                if isinstance(c1,(float,int)) and isinstance(c2,(float,int)):
                    if abs(c1-c2)>self.tol:
                        return False
                else:
                    return False
        return True 
        #return self.values==other.values
        #return all(c1==c2 for c1,c2 in zip(self.values,other.values))
    def __iter__(self):
        return iter(self.values)

    def deriv(self,dervar):
        newcomps=[]
        for comp in self.values:
            if isinstance(comp,(int,float)):
                newcomps.append(0)
            else:
                newcomps.append(comp.deriv(dervar))
        #print(newcomps)
        return Vector(*newcomps)
    def grad(self):
        #print(self,"entering grad(), values:",self.values)
        vars=list(sorted(self.get_vars()))
        
        if vars:
            newcomps=[]
            for comp in self.values:
                if isinstance(comp,(int,float)):
                    newcomps.append(Vector(*(0 for _ in self.values)))
                else:
                    newcomps.append(Vector( *(comp.deriv(vars[ind]) if ind<len(vars) else Monom(0)  for ind,_ in enumerate(self.values))))
                    
            return Vector(*newcomps)
        else:
            return Vector(*(Vector(*(0 for _ in self)) for _ in self))

    def divergence(self):
        vars=list(sorted(self.get_vars()))
        return sum(comp.deriv(var) if isinstance(comp,GeneralObject) else Monom(0) for comp,var in zip(self.values,vars))
    def curl(self):
        if self.dim==2:
            vec=Vector(*self.values,0)
        else:
            vec=self
        
        if vec.dim!=3:
            raise ValueError("The curl can only be taken in 3 dimensions!")
        (_,xy,xz),(yx,_,yz),(zx,zy,_)=( comp.values for  comp in vec.grad())

        return Vector(zy-yz,xz-zx,yx-xy)

    def laplace(self):
        vars=list(sorted(self.get_vars()))
        return Vector(*(sum(comp.deriv(var).deriv(var) if isinstance(comp,GeneralObject) else Monom(0) for var in vars) for comp in self.attached) )
    def simplify(self):
        return self.__class__(*(comp.simplify() if isinstance(comp,GeneralObject) else comp 
                                for comp in self.values))
    @property                        
    def real(self):
        return Vector(*(comp.real if isinstance(comp,(complex,Complex)) else comp for comp in self))
    @property
    def imag(self):
        return Vector(*(comp.imag if isinstance(comp,(complex,Complex)) else 0 for comp in self))
    def flat(self):
        res=[]
        for comp in self:
            if isinstance(comp,Vector):
                for subcomb in comp.flat():
                    res.append(subcomb)
            else:
                res.append(comp)
        return Vector(*res)
    def __abs__(self):
        return self.__class__(*(abs(comp) for comp in self))

class ImplicitVector(Vector):
    def __init__(self,*values):
        #print("Vector got values",values)
        self.attached=self.values=list(it.chain.from_iterable(val if isinstance(val,ImplicitVector) else (val,) for val in values))#list(comp for vec in values for comp in vec)
        #if all(isinstance(comp,ImplicitVector) for comp in values):
        #    self.attached=self.values=list(comp for vec in values for comp in vec)# sum(values,[])
        #else:
        #    self.attached=self.values=list(values)
        #self.attached=self.values=list(values)
        self.dim=len(self.values)
        #self.tol=1.e-10

    def __repr__(self):
        return f"ImplicitVector("+",".join(repr(value) for value in self.values)+")"
    pass
    def __str__(self):
        #print("vvals:",self.values)
        #if len(self.values)==1:
        #    return f"IV({self.values[0]})"
        #else:
            return f" , ".join(truncate_number(value) if isinstance(value,float) else "".join(map(str,value))
        if isinstance(value,tuple) else str(value) for value in self.values)
class VectorGenerator(Func):
    @property
    def name(self):
        return "Vector Generator"
    @property
    def evaluables(self):
        return Vector,int,float
    @property
    def numarg(self):
        return 2
    @property
    def func(self):
        def make_vector(larg,rarg):
            if isinstance(larg,(int,float)):
                if isinstance(rarg,Vector):
                    rest=list(rarg[2:])
                    rarg=rarg[1]
                else:
                    rest=[]
                if  isinstance(rarg,(int,float)):
                    dt=np.int if (type(larg)==type(rarg)==int) else np.float
                    if rarg>=larg:
                        dt=int if (isinstance(larg,int) and isinstance(rarg,int)) else float
                        
                        return ImplicitVector(*it.chain((dt(el) for el in np.arange(larg,rarg+1,dtype=dt)),rest))
                        #else:
                        #    return Vector(*np.concatenate((np.arange(larg,rarg+1),rest)))
                    else:
                        #if type(larg)==type(rarg)==int:
                        return ImplicitVector(*it.chain((dt(el) for el in np.arange(larg,rarg-1,-1,dtype=dt)),rest))
                        #return Vector(*np.concatenate((np.arange(larg,rarg-1),rest)))
                else:
                       return VectorGenerator(larg,rarg)
            else:
                la=len(larg.attached)
                if la==0:
                    return VectorGenerator(0,rarg).eval()
                elif la==1:     
                    return VectorGenerator(larg[1],rarg).eval()
                else:
                    lest=list(larg[:-2])
                    l1,l2=larg[la-1],larg[la]
                    dif=l2-l1
                    if isinstance(rarg,Vector):
                        rest=list(rarg[2:])
                        rarg=rarg[1]
                    else:
                        rest=[]
                    if all(isinstance(el,(int,float)) for el in (l1,l2,rarg)):
                        dt=np.int if all(isinstance(el,int) for el in (l1,l2,rarg)) else np.float
                        #direction=1 if dif>0 else (-1 if dif<0 else 0)
                        #return ImplicitVector(*it.chain(lest,(dt(el) for el in np.arange(l1,rarg+dif,dif,dtype=dt)),rest))
                        
                        return ImplicitVector(*it.chain(lest,(dt(el) for el in np.arange(l1,rarg+dif/1.e9,dif,dtype=dt)),rest))
                    else:
                        return VectorGenerator(larg,rarg)

        return make_vector
    @property
    def inverse(self):
        return type(None)
    def __str__(self):
        return f"{self.attached[0]}...{self.attached[1]}"

    
    

class Abs(Func):
    @property
    def name(self):
        return "abs"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return int,float,Complex,complex
    @property
    def inverse(self):
        return type(None)
    @property
    def func(self):
        def compeval(z):
            if isinstance(z,(Complex,complex)):
                if isinstance(z.real,(int,float)) and isinstance(z.imag,(int,float)):
                    return abs(z.real+1j*z.imag).real
                else:
                    return Abs(z)
            else:    
                return abs(z)
        return compeval
    @property
    def numpyfunc(self):
        return np.abs
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Sign(self.attached[0])*(self.attached[0]).deriv(dervar)

class Angle(Func):
    @property
    def name(self):
        return "angle"
    @property
    def numarg(self):
        return 2
    @property
    def evaluables(self):
        return int,float,Complex
    @property
    def inverse(self):
        return type(None)

   
    @property
    def func(self):
        def angle(rreal,rimag):
            z=Complex(rreal,rimag)
            if isinstance(z.real,(int,float)) and isinstance(z.imag,(int,float)):
                return cmath.phase(z.real+1j*z.imag).real
            else:
                return Angle(rreal,rimag)
        return angle
        #return self.complexification()
    @property
    def numpyfunc(self):
        def angle(real,imag):
            return np.angle(real+1j*imag)
        return angle
    @property
    def ofunc(self):
        def angle(rreal,rimag):
            num=Complex(rreal,rimag)
            real,imag=num.real,num.imag
            
                
            conv=1 #180/math.pi
            if real>0:
                if imag>=0:
                    return math.atan(imag/real)*conv
                else:
                    return (2*math.pi+math.atan(imag/real))*conv
            elif real<0:
                return (math.atan(imag/real)+math.pi)*conv
            else:
                if imag>=0:
                    return math.pi/2*conv
                else:
                    return 3*math.pi/2*conv
        return angle

class Conj(Func):
    @property
    def name(self):
        return "conj"
    @property
    def numarg(self):
        return 1
    @property
    def evaluables(self):
        return Complex,complex
    @property 
    def inverse(self):
        return Conj
    @property
    def func(self):
        def func(z):
            return Complex(z.real,-z.imag)
        return func
    @property
    def numpyfunc(self):
        return np.conj
    def simplify(self):
        arg=self.attached[0]
        if isinstance(arg,self.inverse):
            return arg.attached[0]
        return self
    def deriv(self,dervar):
        if not isinstance(self.attached[0],GeneralObject):
            return 0
        else:
            return Conj(self.attached[0].deriv(dervar))

class Complex(GeneralObject):

    @property
    def tol(self):
        return 1.e-10
    @property
    def name(self):
        return "Complex"
    @property
    def evaluables(self):
        return int,float,complex,Vector
    @property
    def is_evaluable(self):
        return all(comp.is_evaluable if isinstance(comp,GeneralObject) else comp in self.evaluables
                   for comp in (self.real,self.imag))
    def __init__(self,real,imag=None):
        if isinstance(real,(Complex,complex)):
            self.real=real.real
            self.imag=real.imag
            if isinstance(imag,(Complex,complex)):

                self.real-=imag.imag
                self.imag+=imag.real        
            else:
                try:
                    self.imag+=imag
                except TypeError:
                    pass
        else:
            self.real=real
            if isinstance(imag,Complex):
                self.real-=imag.imag
                self.imag=imag.real
            else:
                if imag is None:
                    self.imag=0
                else:
                    self.imag=imag
        self.attached=self.values=[self.real,self.imag]
        self.dim=2 #len(self.values)
    def eval(self,*args,**at):
        try:
            rres=self.real.eval(**at)
        except (AttributeError,TypeError):
            rres=self.real
        try:
            ires=self.imag.eval(**at)
        except (AttributeError,TypeError):
            ires=self.imag
        #return rres+Complex(0,1)*ires
        #print(rres,ires)
        if isinstance(rres,Vector):
            if isinstance(ires,Vector):
                
                return Vector(*(Complex(c1,c2) for c1,c2 in it.zip_longest(rres,ires,fillvalue=0)))
            else:
                return Vector(*(Complex(c1,ires) for c1 in rres))
        elif isinstance(ires,Vector):
            return Vector(*(Complex(rres,c2) for c2 in ires))

        return Complex(rres,ires)
    @property
    def numpyfunc(self):
        def func(x,y):
            return x+1j*y
        return func
    @ignore("Equation")
    def __add__(self,other):
        #print(type(other))
        if isinstance(other,Complex):
            return Complex(self.real+other.real,self.imag+other.imag)
        elif isinstance(other,Vector):
            return Vector(*(self+comp for comp in other.values))            
        else:
            return Complex(self.real+other,self.imag)
    __radd__=__add__
    def __neg__(self):
        return Complex(-self.real,-self.imag)
    @ignore("Equation")
    def __sub__(self,other):
        return self+(-other)
    def __rsub__(self,other):
        return (-self)+other
    @ignore("Equation")
    def __mul__(self,other):
        if isinstance(other,Complex):
          #  #print("compmul",
          #  " | ".join(str(el) for el in (self.real,
           # self.imag,repr(self.real),repr(other.real),repr(self.real*other.real),type(self.real*other.real))))
            return Complex(self.real*other.real-self.imag*other.imag,self.real*other.imag+self.imag*other.real)
        elif isinstance(other,Vector):
            return other*self
        else: 
            #print("mul",repr(self),repr(other),Complex(self.real*other,self.imag*other))
            #print(repr(self.real*other),repr(self.imag*other))
            return Complex(self.real*other,self.imag*other)
    __rmul__=__mul__
    def dot(self,other):
        return self.real*other.real+self.imag*other.imag
    @ignore("Equation")
    def __truediv__(self,other):
        if isinstance(other,Complex):
            return self/other.dot(other)*other.conj()
        elif isinstance(other,Vector):
            return Vector(self/comp for comp in other.values)
        else:
            return Complex(self.real/other,self.imag/other)
    def __rtruediv__(self,other):
        return other/self.dot(self)*self.conj()
    def norm(self):
        return Abs(self).eval()#Sqrt(s elf.real**2+self.imag**2).eval()
    def __abs__(self):
        return self.norm()
    def conj(self):
        return Complex(self.real,-self.imag)
    
    def angle(self):
       # res=Angle(self.real,self.imag).eval()
        return Angle(self.real,self.imag).eval()
    def __pow__(self,exponent):
        ang=self.angle()#*math.pi/180
        #print(ang)
        if isinstance(exponent,int):
            if exponent>=0:
                return ft.reduce(mul,(self for _ in range(exponent)),1)
        return self.norm()**exponent*Complex((Cos(exponent*ang)).eval(),Sin(exponent*ang).eval())
    def __rpow__(self,base):
        if isinstance(base,(int,float)):
            x,y=self.real,self.imag
            return (base**x*(Complex(Cos(y*Ln(Complex(base))),Sin(y*Ln(Complex(base)))))).eval()

    def __str__(self):
        
        real,imag=self.regularized().values
        lpar,rpar=["",""]
        if self==0:
            return "0"
        if real==0:
            if imag==0:
                return "0"
            else:
                rstr=""
        else:
            if isinstance(real,(int,float)):
                rstr=truncate_number(real)
                
            else:
                rstr=str(real)  
            if imag==0:
                return rstr
            try:
                    if imag>0:
                        rstr+=" + "  
                    else:
                        rstr+=" - "  
                        imag*=-1
            except TypeError:
                if isinstance(imag,Monom):
                    if imag.coef>0:
                        rstr+=" + "  
                    else:
                        rstr+=" - "  
                        imag=-imag
                        #imag.coef*=-1
                else:
                    if imag==0:
                        rstr+=""
                    else:
                        rstr+=" + "
        
        if imag==1:
            ipref=""#str(imag)
        elif imag==-1:
            ipref="-"
        else:
            if isinstance(imag,(int,float)):
                ipref=truncate_number(imag)
                #if imag<0:
                #    lpar,rpar="()"
            else:
                ipref=str(imag)
                if  isinstance(imag,(Monom,Sqrt)):
                    return rstr+ipref+"i"
                else:
                    return rstr+"i"+"("+ipref+")"

        return rstr+lpar+ipref+rpar+"i"
    def isreal(self):
        if isinstance(self.imag,(int,float)):
            return abs(self.imag)<self.tol
        else:
            return self.imag==0
    def isimag(self):
        if isinstance(self.real,(int,float)):
            return abs(self.real)<self.tol
        else:
            return self.real==0
    
    def regularized(self):
        real,imag=self
        if isinstance(self.real,(float)):
            if abs(self.real)>self.tol:
                real=self.real
            else:
                real=0
        if isinstance(self.imag,(float)):
            if abs(self.imag)>self.tol:
                imag=self.imag
            else:
                imag=0
        return Complex(real,imag)
        
    def __eq__(self,other):
        if isinstance(other,Complex):
            
            for c1,c2 in zip( (self.real,self.imag),(other.real,other.imag)):
                if c1!=c2:
                    if isinstance(c1,(float,int)) and isinstance(c2,(float,int)):
                        if abs(c1-c2)>self.tol:
                            return False
                    else:
                        return False
            return True 
        else:
            return self.real==other and self.isreal()
    def __float__(self):
        return float(self.real)
        
    def __int__(self):
        return int(self.real)
    def __complex__(self):
        return complex(self.real,self.imag)
    def __repr__(self):
        return f"Complex({repr(self.real)},{repr(self.imag)})"
    def __hash__(self):
        return hash(repr(self))
    def contains_complex(self):
        return True        
    def deriv(self,dervar):
        #return Complex(*super().deriv(dervar))                       
        rder=0 if not isinstance(self.real,GeneralObject) else self.real.deriv(dervar)
        ider=0 if not isinstance(self.imag,GeneralObject) else self.imag.deriv(dervar)
        return Complex(rder,ider)                       
    def __iter__(self):
        return iter(self.values)               
    def simplify(self):
        return self.__class__(*(comp.simplify() if isinstance(comp,GeneralObject) else comp 
                                for comp in (self.real,self.imag)))                      
    
def test():
    x=Monom(3,{"x":2,"y":3})
    y=Monom(9,{"x":4,"y":12})
    z=Monom(12,{"x":4,"y":12,"z":3})
    #print("Factoring ",x+y+z)
    #print("||".join(str(el) for el in (x+y+z).factor_variables()))
    x=Monom(1,{"x":1})
    y=Monom(1,{"y":1})
    z=Monom(1,{"z":1})
    #p=6*x**2*y**3+9*x**4*y**5+15*x**5*y**4*z**3
    #q=6*x**2*y**3+9*x**4*y**5+15*x**5*y**4*z**3
    #z=p/(2*q)
  
    #print("Factoring ",p)
    #print(" | ".join(str(el) for el in p.factor_variables()))
    #print("Factoring ",q)
    #print(" | ".join(str(el) for el in q.factor_variables()))
    #print("Canceling ",z)
    #z=z.cancel_common_prefactors()
    #print(x/y)
    #print("The result is: ",z)

"""
def factor_variables(self):
    cf=1
    d=None
    d=self.monlist[0].coef
    for mon in self.monlist[1:]:
        d=gcd(d,mon.coef)
    pref_vars={var:min(mon.vars.get(var,0) for mon in self.monlist) for var in self.varset}
    prefactor=Monom(d,pref_vars)
    new_coefs=(mon.coef//d for mon in self.monlist)
    new_vars=({var:mon.vars[var]-pref_vars[var] for var in mon.vars}  for mon in self.monlist)
    return (prefactor,Polynomial(Monom(coef,var) for coef,var in zip(new_coefs,new_vars)))
"""


def gcd(a,b,tol=1.e-10):
    #print("gcd:",a,b)
    #print("entering gcd",a,b)
    if a==0:
        return b
    if b==0:
        return a
    
    if not (int(a)==a and int(b)==b):
        if abs(a%b)<tol:
            #print(a%b)
            #print("Warning! a or b is not integer! Since b | a, I return  b")
            return b
        elif abs(b%a)<tol:
            #print("Warning! a or b is not integer! Since a | b, I return  a")
            return a
        else:
            #print("Warning! a or b is not integer! Since not b | a nor a|b , I return  1")
            return 1 # max(abs(a),abs(b))/min(abs(a),abs(b))
    #print("leaving gcd",a,b)
    a,b=(a,b) if a>=b else (b,a)
    while True:
        if b==0:
            return a
        else:
            a,b=b,a%b


def lcm(a,b):
    return a*b//gcd(a,b)

def FFT(a,inverse=False):
    N=len(a)
    power=math.log2(N)
    fac=-2j*np.pi
    if inverse:
        fac*=-1
    if int(power)!=power:
        raise ValueError("Vzorků musí být 2^N!")
    def _FFT(a,N):
        if N==1:
            return np.array(a[0])
        else:
           # print("N",N)
            k=np.arange(N//2)
            wNk=np.exp(fac*k/N)
            even=a[::2]#np.array([a[2*i] for i in k ])
            odd=a[1::2]#np.array([a[2*i+1] for i in k ])
         #   print(np.shape(even),np.shape(odd))
            Seven=_FFT(even,N//2)
            Sodd=_FFT(odd,N//2)
            #print("Seven",Seven,"Sodd",Sodd)
            #print("Sshapes",N//2,np.shape(even),np.shape(odd),np.shape(Seven),np.shape(Sodd))
            #input()
            return np.hstack([Seven+wNk*Sodd,Seven-wNk*Sodd])
    if inverse:
        return _FFT(np.array(a),N)
    else:
        return _FFT(np.array(a),N)/N

def large_binom(p,n,k):
    res=1
    d=k    
    if k>n:
        raise ValueError("Can't evalute binom(p,n,k) for n,k integers and n<k.")
    if k<d:
        k=n-k
        d=n-k
    if (k>=d):
        for q in range(n,k,-1):
            
            res*=q/(q-k)*(1-p)
            #print(q,q-k,d,n,res)
            while res>1 and d>=0:
                    res*=p
                    d-=1
        return res*p**(d)
        
#if __name__=="__main__":
#    x=Variable("x")
#    y=Variable("y")
#    z=Variable("z")
    
    #print( (1/Cos(x))/Cos(x))
    """x=Variable("x")
    y=Variable("y")
    #q=x*y
    #s=Sin(x*y)
    #z=Sin(x*y)
    #w=Monom(1,{str(s):(s,1)})
    #a=5*(x+y)*x
    #a=5*z+s
    #b=s+5*z+Cos(x)
    #print(a,b,a.eval(y=1),b.eval(y=1),a.eval(x=1,y=5),b.eval(x=1,y=7))
    #print((x+y)**3,"|||",((x+y)**3).deriv("x"))
    #print((Cos(x*y**2)+y+y**2*Sin(x)).deriv("y"))
    v=(Sin(x)*Cos(x)*x*y).deriv("x")
    #print(v,v.eval(x=1))
    #print(*(repr(vy) for vy in v.monlist))
    #print(Sin(4)*Sin(5**2)**5)
    #print(Sin(x)**(x)) 
    #print(Sin(x),Sin(x).eval(x=5),(x+Sin(x)**x)**2,((x+Sin(x)**x)**2).eval(x=4))
    #print(((Sin(x)/y).eval(x=Sin(x)**2)).eval(y=100,x=3),x/x,(x*y)/(x**2))
    x=Variable("x")
    y=Variable("y")
    #print(x.eval(x=4),x**2,y.eval(y=7))
    #print(x.name,repr(x),str(x))
    #print(x/y,(x**2)/y,(Sin(x)/Cos(y)).eval(x=Cos(y),y=Cos(y)))
    #print(x.eval(y=4),(y*x**2).eval(y=4,x=4),x*y,x*5*Sin(x)**2)
    #print(s,s.eval(x=1,y=2))
    #print(q,q.eval(x=1,y=5))
    #print(s*q,(s*q).eval(x=1,y=2))
    #print(s*q*s,(s*q*s).eval(x=1))
    #print(w*w)
    #print(x)
    #print(x.eval(x=3))
    #print(s)
    #q=Monom(4,{"("+str(s)+")":(s,2)})
    #s=go.Sin(go.Cos(x*q)) 
    #s=go.Sin(go.Sin(x*y))
    #print(y,z)
    
    #print(s)
    #print("...")
    #print(s.eval(),s.eval(x=2),s.eval(y=2),s.eval(x=2,y=2))
    #print(y.eval(x=0))
    #q=Monom(1,{"("+str(x)+")":(x,2)})
    #print(q,q**2)
    #print(q.eval(y=2),(q**2).eval(y=2))
    #s=go.Sin()
    x=Variable("x")
    y=Variable("y")
    s=Sin(x)
    #print(x,x**2,x+3*y,(x+3*y)**2,s,s**2,(s+x)**2,(s+x)**x,((s+x)**x)**3)
    vars=(x,x**2,x+3*y,(x+3*y)**2,s,s**2,(s+x)**2,(s+x)**x,((s+x)**x)**3,x.deriv("x"),s.deriv("x"))
    #print(" | ".join(str(var.eval()) for var in vars))
    #print(" |".join(str(var.eval(x=2)) for var in vars))
    #print(" | ".join(str(var.eval(y=2)) for var in vars))
    #print(" |".join(str(var.eval(x=2,y=2)) for var in vars))
    z=Sin(x)*Cos(y)**2
    d=(z.deriv("y")).deriv("x")
    #print(z,z.deriv("y"),z.deriv("y").deriv("x"),z.deriv("y").deriv("x").eval(y=x,x=y))
    z=Sin(2*x)
    #print(z,z.simplify())
    q=Cos(2*x)
    #print(q,q.simplify())
    t=Tan(2*x)
    #print(t,t**2,t.simplify())
    y=Sin(Acos(x))
    #print(y,y.simplify(),y.eval(x=1))
    for func in (Sin,Cos,Tan,Asin,Acos,Atan,Sqrt):
        #print(func(x).deriv("x"),end="|")
    z=Variable("z")
    #print("\n",y.get_vars(),x.get_vars(),Sin(z+x).get_vars())
    #print((x+z*z+Sin(x)).grad())
    #print((y*y).deriv("x"))
    #print(Vector(x,z).divergence())
    y=Variable("y")
    lap=(1/Sqrt(x**2+y**2+z**2)).laplace()
    #print( (1/Sqrt(x**2+y**2+z**2)).grad().divergence())
    #print( (1/Sqrt(x**2+y**2+z**2)).laplace())
    #print(Vector(x*z,y*z).divergence())
    #print(lap," | ",lap.simplify())
    #print(Sqrt(x)**2,(3*Sqrt(x+y)**2).simplify())
    #print(((x**(x))**(1/x)).simplify())
    #print(Ln(x**2).deriv(x))
    #print(1/x**2*2*x)
    #print((y*y).deriv("y"))
    #print(Sin(x).eval(x=1),1/(1+Sin(x).eval(x=1)))
    #print(z,1/(z+1),z.eval(x=4),(1/z).eval(x=4))"""
    """Angle, old eval:
     def eval(self,*args,**kwargs): 
        if len(self.attached)!=self.numarg:
            raise ValueError(f"Func {self.name}  needs to take exactly {self.numarg} argument(s)!")
        else:
            if Vector not in self.evaluables and type(self.attached[0]) is Vector:
                return Vector(*(self.__class__(comp).eval(*args,**kwargs) for comp in self.attached[0]))
            z=Complex(self.attached[0],self.attached[1])
            results=[obj.eval(*args,**kwargs) if isinstance(obj,GeneralObject) else obj for obj in z]
            if all(isinstance(res,self.evaluables) for res in results):
                return self.func(*(results))
            else:
                return self.__class__(*results)
    """