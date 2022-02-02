import math
from .gen_obj import Polynomial,SPolynomial,Vector,Quotient,truncate_number,Sqrt,Complex,Cos,Sin

def solve_linear(a,b,fractions=True):
    num=-b
    den=a

    if den==0:
        if num==0:    
            return "Always true"
        else:
            return "No solution"
    else: 
        if fractions:
            ##print("fractions")
            #num,den=make_whole(num,den)
            return Quotient(num,den)
        else:
            return num/den if num%den else int(num//den)

def solve_quadratic(a,b,c,fractions=True,format=True,use_complex=True):
    #print("fractoins,format,compl:",fractions,format,use_complex)
    
    if a==0:
        return solve_linear(b,c,fractions=fractions)
    bh=(b)/2
    D=bh**2-a*c
    iD=int(D)
    if D==iD:
        D=iD
    #print("a-D",a,bh,c,D)
    imag=False
    #print(D,type(D))
    if isinstance(D,Complex):
        #print("D is complex")
        sqrtD1,sqrtD2=roots(D,2)
        if fractions:
            return Vector(Quotient(-bh+sqrtD1,a),Quotient(-bh+sqrtD2,a))
        else:
            return Vector(-bh+sqrtD1,-bh+sqrtD2)/a
    if D<0:
        if not use_complex:
            return "No solution"
        else:
            imag=True
    elif D==0:
        if fractions:
            return Vector(Quotient(-bh,a),Quotient(-bh,a))
            
        else:
           return Vector(-bh/(a),-bh/a)
    
    if imag:
        D=-D
    sqrtD= math.sqrt(D)    
    if fractions:
        
        #print("fractions")
        if int(sqrtD)==sqrtD or int(1/sqrtD)==1/sqrtD:
            #print("got here",sqrtD,type(sqrtD))
            if imag:
                sqrtD=Complex(0,sqrtD)
            num1,num2,den=-bh+sqrtD,-bh-sqrtD,a
            print("here",num1,num2,den,"|",repr(num1),repr(num2),repr(den))
            return Vector(Quotient(num1,den),Quotient(num2,den))  #bacha!
        vertex=Quotient(-bh,a)
        
        if format:
            if imag:
                return Vector(-bh+Complex(0,Sqrt(D)),-bh+Complex(0,-Sqrt(D)))
            else:
                return Vector(-bh+Sqrt(D),-bh-Sqrt(D))
            D=truncate_number(D)
            
            #num,den=make_whole(-bh,a)
            #vertex=Quotient(num,den)
            if imag:
                prefix="i"
            else:
                prefix=""
            ##print(vertex,repr(vertex),vertex.p.isconst,vertex.q.isconst)
            if vertex==0: 
                        vertex=""
            if a<0:
                a=-a
            if a==1:
                return (Vector((vertex,f"+{prefix}√({D})"),(vertex,f"-{prefix}√({D})")))
            else:
                return (Vector((vertex,f"+{prefix}√({D})/{a}"),(vertex,f"-{prefix}√({D})/{a}")))
        else:
            gensqrt=Sqrt(D)
            if imag:
                gensqrt*=Complex(0,1)
            return Vector( Quotient( (-bh+gensqrt).eval(),a),Quotient((-bh-gensqrt).eval(),a))
    
    else:
        if imag:
            sqrtD*=Complex(0,1)
        return (Vector(sqrtD,-sqrtD)-bh)*(1/a)

def solve_cubic(a,b,c,d,fractions=True,use_complex=True):
    tol=1.e-10  
    ##print(a,b,c,d)
    if a==0:
        return solve_quadratic(b,c,d,fractions=fractions,use_complex=use_complex,format=False)
    
    if d==0:
        qsol=solve_quadratic(a,b,c,fractions=fractions,format=False,use_complex=use_complex)
        ##print(qsol)
        if qsol=="No solution":
            return Vector(0)
        else:
            return Vector(0,*qsol)
    else:
        ##print(a,b,c,d)
        b,c,d=b/a,c/a,d/a
        ##print(b,c,d)
        if b!=0:
            p,q=depress(b,c,d)
        else:
            p,q=c,d
        ##print("p","q",p,q)
        if p==q==0:
            t1,t2,t3=0,0,0
        elif p==0:
            t1,t2,t3=regroots(-q,3)
        elif q==0:
            t1=0
            t2,t3=regroots(-p,2)
        else:
            D=(q/2)**2+(p/3)**3
            ##print(a,b,c,d,q,p,"D:",D)
         #   #print(D)
            if abs(D)<tol:
               # #print(-q,-q/2)
                u=realroots(-q/2,3)
                #v=realroots(q/2,3)
                ##print("D=0",u,b/3)
                if u:   
                    t1=2*u[0]
                   # t1=u[0]-v[0]
                    t2,t3=-u[0],-u[0]
                else:
                    t1,t2,t3=roots(-q/2,3)
            else:
                
                sqrtD=regroots(D,2)[0]
               

                u1,u2,u3=roots(-q/2+sqrtD,3)
                #v3,v1,v2=roots(q/2+sqrtD,3)

                #t1,t2,t3=(u1-v1,u2-v2,u3-v3)
                t1,t2,t3=tuple((u-p/(3*u)).regularized() for u in (u1,u2,u3))
                ##print(" | ".join(str(el) for el in (u1,u2,u3,v1,v2,v3,t1,t2,t3,D,sqrtD)))
            
        x1,x2,x3=( (t-b/3) for t in (t1,t2,t3))
        ##print(x1,x2,x3)
        if use_complex:
            return Vector(*(x.regularized() if isinstance(x,Complex) else x for x in  (x1,x2,x3)))
        else: 
            return Vector(*(x for x in (x1,x2,x3) if not isinstance(x,Complex) or x.isreal() ))
            #t**3+pt+q=0

def solve_quartic(a,b,c,d,e,fractions=True,use_complex=True):
    
    if a==0:
        return solve_cubic(b,c,d,e,fractions=fractions,use_complex=use_complex)
    if e==0:
        return Vector(0,*solve_cubic(a,b,c,d,use_complex=use_complex,fractions=fractions))
    
    b,c,d,e=(num/a for num in (b,c,d,e))
    p,q,r=depress_quartic(b,c,d,e)
    #print("bcde:",b,c,d,e)
    #print("pqr:",p,q,r)
    if q==0:
        y1,y2=solve_quadratic(1,p,r,fractions=False,use_complex=True,format=False)
        x1,x2=(x.regularized() for x in roots(y1,2))
        x3,x4=(x.regularized() for x in roots(y2,2))
        return Vector(*(sol-b/4 for sol in (x1,x2,x3,x4)))

    res_sols=solve_cubic(*resolvent_cubic(p,q,r),fractions=False,use_complex=True)                                                                          
    for sol in res_sols:
        if sol!=0:
            m=sol
            if  (not isinstance(sol,Complex) or sol.isreal):
                m=float(sol)
                break
    
    sm1,sm2=roots(2*m,2)
    y1,y2=solve_quadratic(1,sm1,p/2+m-q/(2*sm1))
    y3,y4=solve_quadratic(1,sm2,p/2+m-q/(2*sm2))
    sols=(y-b/4 for y in (y1,y2,y3,y4))
    return Vector(*(x.regularized() if isinstance(x,Complex) else x for x in sols))
        

def depress(b,c,d):
    return (3*c-b**2)/3,(2*b**3-9*b*c)/27+d
def depress_quartic(b,c,d,e):
    return -(3*b**2)/8+c,(b**3-4*b*c)/8+d,(-3*b**4-64*b*d+16*b**2*c)/256+e
def resolvent_cubic(p,q,r):
    return 8,8*p,2*p**2-8*r,-q**2
def roots(num,n):
   z=Complex(num)
   norm,ang=abs(z)**(1/n),z.angle()/n 
   #print(norm,ang)
   res=[norm*(Cos(arg).eval()+Complex(0,1)*Sin(arg).eval()) for arg in (ang+i/n*2*math.pi for i in range(n))]
   #print(res)
   return tuple(res)

def regroots(num,n):
    return tuple(root.real if root.isreal() else root for root in roots(num,n))

def realroots(num,n):
    return tuple(root.real for root in roots(num,n) if root.isreal())