from . import gen_obj
from .gen_obj import Complex,Variable,Vector
import random
import math,cmath
def roots_key(root):
    return abs(root)+cmath.phase(root)

def absolve(pol):
    try:
        
        var=next(iter(pol.varset))
        varname=next(iter(pol.get_vars()))
        coefs=pol.get_coefs(var)
        lbound,ubound=get_bounds(coefs)
      #  print(lbound,ubound)
        if ubound==0:
            return Vector(*(0 for _ in coefs[:-1]))
        boundlog=round(math.log10(ubound))
        niters=20 +(boundlog if boundlog<3 else boundlog*10)
    
        pderiv=pol.deriv(varname)
        ratio=pol/pderiv
        for _ in range(100):
            zks=gen_compl(ubound,pol.order)
         #   print(pol,repr(pol),"||",varname)    
            #print(repr(pol),"||",repr(Variable(var)))
            #print("Try",zks,var,coefs,pderiv)
            for _ in range(niters):
                qks=[ratio.eval(**{varname:zk}) for zk in zks]
                zks=[qk.real+1j*qk.imag for qk in zks]
                qks=[qk.real+1j*qk.imag for qk in qks]
                denoms=[1-qk*sum(1/(zks[k]-zj) for j,zj in enumerate(zks) if j!=k ) for k,qk in enumerate(qks)]
                #wks=[ (qk/(denom)) for zk,qk,denom in zip(zks,qks,denoms)]
                zks=[ zk-(qk/(denom)) for zk,qk,denom in zip(zks,qks,denoms)]
                #print("wk:","|".join(str(w) for w in wks))
                #zks=[zk-wk for zk,wk in zip(zks,wks)]
                #zks=[Complex(zk) for zk,wk in zip(zks,wks)]

                if any(abs(zk)>ubound for zk in zks):
                    break
            #num=pol.leading_term()*mul(z-Complex(zk))
            #print(zks)
            if all(abs(zk)<ubound for zk in zks):
                break
        
        return Vector(*sorted((Complex(zk).regularized() for zk in zks),key=roots_key))
    except ZeroDivisionError:
        #print("zdiver")
        raise ValueError("Division by zero in the approximate solving algorithm, try again")
    except KeyError:
        if len(pol.varset)==0:
            return "Never True" if pol!=0 else "Always True"
        else:
            raise NotImplementedError("Can only solve for functions of one variable")
        
def get_bounds(coefs):
    an=coefs[-1]
    ubound=1+max(abs(a_i/an) for a_i in coefs[:-1] )
    #print(coefs)
    k=len(coefs[:-1])
    #print(k)
    expons=[((k-i)) for i,a_i in enumerate(coefs[:-1])]
    #print(expons,an)
    
    fujiwara_bound=2*max(abs(a_i/an)**(1/(expons[i])) for i,a_i in enumerate([coefs[0]/2]+coefs[1:-1]) )
    #print("uf",ubound,fujiwara_bound)
    ubound=min(ubound,fujiwara_bound)

    coefs=[coef for coef in coefs[::-1] if coef!=0] #reverse
    an=coefs[-1]
    try:
        lbound=1/(1+max(abs(a_i/an) for a_i in coefs[:-1] ))
    except ValueError:
        lbound=0
    return lbound,ubound
def gen_compl(bound,n):
    res=[]
    for _ in range(n):
        norm=random.random()*bound
        rind=random.random()*2*math.pi
        res.append(complex(norm*math.cos(rind),norm*math.sin(rind)))
    return res

