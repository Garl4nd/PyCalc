from . import gen_obj
from .gen_obj import * 

class DifEq(GeneralObject):
    @property
    def name(self):
        return "DifEq"
    @property
    def evaluables(self):
        return int,float
    @property
    def is_evaluable(self):
        return True

    def __init__(self,*args,t0=0.0,t1=1.0,n=100,is_map=False,bifurc=False):
        self.is_map=is_map
        try:
            #self.args=Vector(*args)
            
            self.args=Vector(*args)
            self.res=None
            self.bifurc=bifurc
            
            args=iter(args)
            self.funcs=[];self.init_points=[]
            try:
                while True:
                    self.funcs.append(next(args))
                    self.init_points.append(next(args))
            except StopIteration:
                if len(self.init_points)<len(self.funcs):
                    raise ValueError("Must supply initial conditions for all functions!")
            
        except TypeError:
            funcs,init_points=self.args
            self.funcs=[funcs]
            self.init_points=[init_points]
        self.t0=t0
        self.t1=t1
        self.n=n
        self.data=list(zip(self.funcs,self.init_points))
        self.attached=self.args

    def eval(self,*args,**kwargs):
        dargs=(arg.eval(*args,**kwargs) if isinstance(arg,GeneralObject) else arg for arg in self.args)
        return DifEq(*dargs,t0=self.t0,t1=self.t1,n=self.n,is_map=self.is_map,bifurc=self.bifurc)
    def calculate(self,t0=0,t1=1,n=100):
        return self.pairwise_dif_eq(self.funcs,self.init_points,t0=t0,t1=t1,n=n)

    def dif_eq(self,func,init_point,t0=0,t1=1,n=100):
        from scipy.integrate import odeint
        try:
            if n<0 or int(n)!=n:
                raise ValueError("n (tden)  must be a non-negative natural number!")
        except TypeError:
            raise ValueError("n (tden) must be a non-negative natural number!")
        if not isinstance(func,GeneralObject):
            func=Monom(func)
        if not isinstance(func,Vector):
            try:
                if isinstance(func,Complex):
                    func=Vector(func)
                else:
                    func=Vector(*func)
            except TypeError:
                func=Vector(func)
        if not isinstance(init_point,Vector):
            
            try:
                if isinstance(init_point,Complex):
                    init_point=Vector(init_point)
                else:
                    init_point=Vector(*init_point)
            except TypeError:
                init_point=Vector(init_point)
        if func.dim!=init_point.dim:
            raise ValueError("The dimensions of the mapping function and of the initial point must be equal!")
        
        if self.is_map:
            p=init_point
            res=[]
            if not self.bifurc:
                res.append(p)
            vars=sorted(func.get_vars())
            for ind in range(n):
                p=func.eval(**{var:comp for var,comp in zip(vars,p)})
                #raise ValueError(self.bifurc)
                if (not self.bifurc) or ind>4*n/5:
                    res.append(p)
                #print(p)
            return Vector(*res),vars
        else:
            tarr=np.linspace(t0,t1,n)
            adjvars=func.get_vars().difference({"t"})
            ordvars={var for var in adjvars if not (len(var)>1 and var.startswith("d")) and var.islower() }     
            dvars=adjvars.difference(ordvars)
            for var in dvars:
                if var.startswith("d") and var[1:] not in ordvars:
                    ordvars.add(var[1:])
                elif var.isupper() and var.lower() not in ordvars:
                    ordvars.add(var.lower())

            
            vars=sorted(ordvars)+sorted(dvars)
            candlist=it.chain(("x","y","z"),(chr(i) for i in range(ord("a"),ord("w"))))
            while len(vars)<len(init_point):
                if "y" in ordvars and "x" not in ordvars:
                    vars=["x"]+vars
                s=next(candlist) # Nebo se musí zavolat iter?
                while s in vars:
                    s=next(candlist)
                vars.append(s)
                    
                
                #ordvars={"y"}
            
            def f(y,t):
                return [*func.eval(**{var:comp for var,comp in zip(["t"]+vars,[t]+list(y))})]
            #print("ip",init_point)
            return np.c_[tarr,odeint(f,init_point,tarr,args=tuple())],["t"]+vars
        #return Vector(*res)
    
    def pairwise_dif_eq(self,funcs,init_points,t0=0,t1=1,n=100):
        #args=list(args)
        #print(args)
        """
        try:
            funcs,init_points=zip(*args)
        except TypeError:
            funcs,init_points=args
            funcs=[funcs]
            init_points=[init_points]
        """
        names=[f"f = ${func}$, init: {init}" for func,init in zip(funcs,init_points)]
        res_array=[]
        
        for func,init_point in zip(funcs,init_points):
            res,var=self.dif_eq(func,init_point,t0=t0,t1=t1,n=n)
            res_array.append(res)
            #res_array.append(self.dif_eq(func,init_point,t0=t0,t1=t1,n=n))
        self.res=res_array[0]
        return res_array,names,var
        
    def __repr__(self):
        return f"ODE object with data: {self.args}, is_map = {self.is_map}"
    def __str__(self):
        if self.res is not None:
            return str(self.res)
        if self.is_map:
            return "System of maps x_n->f(x_n) with: \n\t"+"\n\t".join(f'f= {func}, init point: {init}' for func,init in zip(self.funcs,self.init_points))
        else:
            return "System of ODEs y'=f with: \n\t"+"\n\t".join(f'f= {func}, init point: {init}' for func,init in zip(self.funcs,self.init_points))

class BVP(DifEq):
    @property
    def name(self):
        return "BVP"
    @property
    def evaluables(self):
        return int,float
    @property
    def is_evaluable(self):
        return True

    def __init__(self,*args,t0=0.0,t1=1.0,n=100,is_map=False):
        self.is_map=is_map
        try:
            #self.args=Vector(*args)
            
            self.args=Vector(*args)
            args=iter(args)
            self.funcs=[];self.init_points=[]
            try:
                while True:
                    self.funcs.append(next(args))
                    ip=next(args)
                    self.init_points.append(ip)
                    if any(isinstance(c,Vector) for c in ip):
                        if any(not isinstance(c,Vector) or len(c)!=3 for c in ip):
                            raise ValueError("Either all or no conditions must be specified by a three-component vector!")

            except StopIteration:
                if len(self.init_points)<len(self.funcs):
                    raise ValueError("Must supply initial conditions for all functions!")
            
        except TypeError:
            funcs,init_points=self.args
            self.funcs=[funcs]
            self.init_points=[init_points]
        
        self.t0=t0
        self.t1=t1
        self.n=n
        self.data=list(zip(self.funcs,self.init_points))
        self.attached=self.args

    def eval(self,*args,**kwargs):
        dargs=(arg.eval(*args,**kwargs) if isinstance(arg,GeneralObject) else arg for arg in self.args)
        return BVP(*dargs,t0=self.t0,t1=self.t1,n=self.n,is_map=self.is_map)
    def calculate(self,t0=0,t1=1,n=100):
        return self.pairwise_dif_eq(self.funcs,self.init_points,t0=t0,t1=t1,n=n)

    def dif_eq(self,func,init_point,t0=0,t1=1,n=100):
        from scipy.integrate import solve_bvp
        try:
            if n<0 or int(n)!=n:
                raise ValueError("n (tden)  must be a non-negative natural number!")
        except TypeError:
            raise ValueError("n (tden) must be a non-negative natural number!")
        if not isinstance(func,GeneralObject):
            func=Monom(func)
        if not isinstance(func,Vector):
            try:
                if isinstance(func,Complex):
                    func=Vector(func)
                else:
                    func=Vector(*func)
            except TypeError:
                func=Vector(func)
        if not isinstance(init_point,Vector):
            
            try:
                if isinstance(init_point,Complex):
                    init_point=Vector(init_point)
                else:
                    init_point=Vector(*init_point)
            except TypeError:
                init_point=Vector(init_point)
        #print(func,init_point)
        if func.dim!=init_point.dim:
            raise ValueError("The dimensions of the mapping function and of the initial point must be equal!")
        
        tarr=np.linspace(t0,t1,n)
        adjvars=func.get_vars().difference({"t"})
        ordvars={var for var in adjvars if not (len(var)>1 and var.startswith("d")) and var.islower() }
        dvars=adjvars.difference(ordvars)
        for var in dvars:
            if var.startswith("d") and var[1:] not in ordvars:
                ordvars.add(var[1:])
            elif var.isupper() and var.lower() not in ordvars:
                ordvars.add(var.lower())

        
        vars=sorted(ordvars)+sorted(dvars)
        candlist=it.chain(("x","y","z"),(chr(i) for i in range(ord("a"),ord("w"))))
        
        while len(vars)<len(init_point):
            if "y" in vars and "x" not in vars:
                vars=["x"]+vars
            else:
                s=next(candlist) # Nebo se musí zavolat iter?
                while s in vars:
                    s=next(candlist)
                vars.append(s)
        #print(init_point)

        init_point=list(init_point)
        def f(t,y):
            #return func.numcompose(y,t)
            #return [*func.eval(**{var:comp for var,comp in zip(["t"]+vars,[t]+list(y))})]
            #print(np.shape(t))
            
            return [*func.numcompose(**{var:comp for var,comp in zip(["t"]+vars,[t]+list(y))})]
        def simple_boundary(y0,y1):
            if len(y0)%2==0:
                n=len(y0)//2
                return [y0[i]-init_point[i] for i in range(n)]+[y1[i]-init_point[i+n] for i in range(n)]
            else:
                return [y0[i]-init_point[i] for i in range(len(y0))]
        def exact_boundary(y0,y1):
            y=[y0,y1]
            rl=[]
            
            for ind,pos,value in init_point:
                rl.append(y[pos][ind-1]-value)
            #print(y0,y1,rl)
            return rl
        #print("ip",init_point)
        if isinstance(init_point[0],Vector):
            boundary=exact_boundary 
            starts,ends=[0 for _ in init_point],[tarr[-1  ] for _ in init_point]
            #y_init=[tarr  for _ in init_point]
            
            y_init=[]
            for ind,pos,value in init_point:
                if pos==0:
                    starts[ind-1]=value
                elif pos==1:
                    ends[ind-1]=value
            for start,end in zip(starts,ends):
                y_init.append(np.linspace(start,end,len(tarr)))
            #print(vars)
            #print(vars,y_init)
            
        else:
            boundary= simple_boundary
            y_init=[tarr  for _ in init_point]

        solver=solve_bvp(f,boundary,tarr,y_init,max_nodes=10000)
        if solver.status!=0:
            print(solver.status)
            raise ValueError(f"THe BVP failed to converge! Vars: {vars}")
        #print(np.shape(tarr),solver.sol(tarr))
        #return tarr,solver.sol(tarr),["t"]+vars
        return np.c_[tarr,np.transpose(solver.sol(tarr))],["t"]+vars

        #return Vector(*res)
    
    def pairwise_dif_eq(self,funcs,init_points,t0=0,t1=1,n=100):
        #args=list(args)
        #print(args)
        """
        try:
            funcs,init_points=zip(*args)
        except TypeError:
            funcs,init_points=args
            funcs=[funcs]
            init_points=[init_points]
        """
        names=[f"f = ${func}$, BC: {init}" for func,init in zip(funcs,init_points)]
        res_array=[]
        
        for func,init_point in zip(funcs,init_points):
            res,var=self.dif_eq(func,init_point,t0=t0,t1=t1,n=n)
            res_array.append(res)
            #res_array.append(self.dif_eq(func,init_point,t0=t0,t1=t1,n=n))
        return res_array,names,var
        
    def __repr__(self):
        return f"ODE object with data: {self.args}, is_map = {self.is_map}"
    def __str__(self):
        if self.is_map:
            return "System of maps with equations: \n\t"+"\n\t".join(f'f= {func}, init point: {init}' for func,init in zip(self.funcs,self.init_points))
        else:
            return "System of BVPs y'=f with: \n\t"+"\n\t".join(f'f= {func}, boundary conditions: {init}' for func,init in zip(self.funcs,self.init_points))