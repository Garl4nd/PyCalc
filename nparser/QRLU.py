import numpy as np
from operator import mul
import functools as ft

tol=1.e-10
def LUdecomp(A):
    U=np.array(A)
    nx,_=np.shape(U)
    L=np.eye(nx)
    perms=[]
    for n in range(nx-1):
        if abs(U[n,n])<tol:
        #    print("zero")
            for i in range(n+1,nx):
                if abs(U[i,n])>tol:
                    U[[i,n]]=U[[n,i]]
                    perms.append((n,i))
                    break
                        #print(perms)
                        #input()

                  #  P[i,i],P[n,n]=0,0
                  #  P[i,n],P[n,i]=1,1
                    
            else:
                li=0
                L[n,n+1:]=li
                continue
    #print("permuted U:",U)
    for n in range(nx-1):
        li=U[n+1:,n]/U[n,n]
        #for i in range(n+1,nx):
        #    U[i,:]=U[i,:]-li[i-n-1]*U[n,:]
        U[n+1:,:]=U[n+1:,:]-np.outer(li,U[n,:])
        L[n+1:,n]=li
    
    return L,U,perms

def LUPsolve(A,Borig):
    B=np.array(Borig)
    L,U,perms=LUdecomp(A)
    lx,*ly=np.shape(B)
    
    if not ly:
        B=np.reshape(B,(lx,1))
        ly=1
    else:
        ly=ly[0]
    x=np.empty((lx,ly))
    y=np.empty((lx,ly))
  #  print("perms",perms)
  #  print("B",B)
    for perm in perms:
        B[[perm[0],perm[1]]]=B[[perm[1],perm[0]]]
    #print("B",B)
    for i in range(lx):
        for k in range(i):
            B[i,:]=B[i,:]-L[i,k]*y[k,:]
        y[i,:]=B[i,:]/L[i,i]
    #print("U:",U,"y:",y)
    for i in reversed(range(lx)):
        for k in range(i+1,lx):
            y[i,:]=y[i,:]-U[i,k]*x[k,:]
        x[i,:]=y[i,:]/U[i,i]
    return x
def LUPinverse(A):
    return LUPsolve(A,np.eye(len(A)))

def QRdecomp(A):
    vs=[vec for vec in np.transpose(A)]
    nx,ny=np.shape(A)
    Q=np.zeros((nx,nx))
    R=np.zeros((nx,ny))
    es=[]
    #print(nx,ny,len(vs),A)
    r=0
    for i in range(ny):
        v=vs[i]
        magv=np.sqrt(abs(np.dot(v,v)))
        if abs(magv)>tol:
            e=v/magv
            es.append(e)
            
            R[r,i]=magv
            for k in range(i,ny):
                R[r,k]=np.dot(vs[k],e)
                vs[k]=vs[k]-R[r,k]*e
            r+=1
            if r==nx:
                break
            
        
    

    
    r=len(es)
    for i in range(r,nx):
        while True:
            v=np.random.random(nx)
            u=v
            #R[:i,i]=0
            for ind,e in enumerate(es):
             
                u=u-np.dot(u,e)*e
            magu=np.sqrt(abs(np.dot(u,v)))
            if magu>=tol:
                es.append(u/magu)
                break
    for i,e in enumerate(es):
        Q[:,i]=np.where(abs(e)<tol,0,e)
    return Q,R,r   

def QRdecomp2(A):
    vs=[vec for vec in np.transpose(A)]
    nx,ny=np.shape(A)
    Q=np.zeros((nx,nx))
    R=np.zeros((nx,ny))
    es=[]
    for i in range(ny):
        u=vs[i]
        #R[:i,i]=[np.dot(u,e) for e in es] 
        for ind,e in enumerate(es):
            R[ind,i]=np.dot(u,e)
            #print("u",u)
            u=u-R[ind,i]*e
        magu=np.sqrt(abs(np.dot(u,vs[i])))
        if magu<tol:
            #print(f"nullifying {i}-th column")
            magu=0
        else:
            es.append(u/magu)
            R[i,i]=magu
    
    
    r=len(es)
    for i in range(r,nx):
       
        while True:
            v=np.random.random(nx)
            u=v
            #R[:i,i]=0
            for ind,e in enumerate(es):
             
                u=u-np.dot(u,e)*e
            magu=np.sqrt(abs(np.dot(u,v)))
            if magu>=tol:
                es.append(u/magu)
                break
  
    for i,e in enumerate(es):
        Q[:,i]=np.where(abs(e)<tol,0,e)
    return Q,R,r   
def proj(v,es):
    #print(np.shape(v),np.shape(es))
    res=sum(np.dot(v,e)/np.dot(e,e)*e for e in es)
    return res

def QRsolve(A,Bc):
    B=np.array(Bc,dtype=float)
    exact=[]
    Q,R,r=QRdecomp(np.transpose(A)) #r = rank
    
    L=np.transpose(R)
    bx,*by=np.shape(B)
    ax,ay=np.shape(A)
    if not by:
        B=np.reshape(B,(bx,1))
        by=1
    else:
        by=by[0]
    y=np.zeros((ay,by))
    ker=Q[:,r:]
    Q2,_,r2=QRdecomp(A)
    if ax==ay==r:
        exact=[True for _ in range(by)]
    else:
        for i,b in enumerate(np.transpose(B)):
            bproj=proj(b,[e for e in np.transpose(Q2[:,:r2])])
            if np.dot(b-bproj,b-bproj)<tol:
                exact.append(True)
            else:
                exact.append(False)
                B[:,i]=bproj
    for i in range(r):  # pivoting
        if abs(L[i,i])<tol:
            for k in range(i+1,ay+1):
                if abs(L[k,i])>tol:
                    L[[i,k]]=L[[k,i]]
                    B[[i,k]]=B[[k,i]]
                    break
    for i in range(r):
            for k in range(i):
                B[i,:]=B[i,:]-L[i,k]*y[k,:]    
            
            y[i,:]=B[i,:]/L[i,i]

    return np.matmul(Q,y),exact,ker

def eigenproblem(Ao): #symmetrizes the problem first, so as to not run into any trouble
    etol=np.sqrt(tol)
    A=np.array(Ao+np.transpose(Ao))/2
    i=0
    counter=0
    last_m=-2*tol
    while True:
        Q,R,_=QRdecomp(A)
        A=np.matmul(R,Q)
        D=diag(A)
        i+=1
        if counter==1000:
            return [D[i,i] for i,_ in enumerate(D)] ,Q
        else:
            m=max(gershgorin(A))

            if m<etol:
                
                return [D[i,i] for i,_ in enumerate(D)] ,Q
            else:
                if abs(m-last_m)<tol:
                    counter+=1
                else:
                    counter=0
                last_m=m

def gershgorin(A):
    r=[]
    for i,row in enumerate(A):
        r.append(sum(abs(el) for k,el in enumerate(row) if k!=i  ))
    return r
def diag(A):
    D=np.zeros(np.shape(A))
    for i,_ in enumerate(A):
        D[i,i]=A[i,i]
    return D
def QRinv(A):
    return QRsolve(A,np.eye(len(A)))
def QRdet(A):
    _,R=QRdecomp(A)
    try:
        return ft.reduce(mul,[R[i,i] for i in range(len(np.transpose(R)))])
    except IndexError:
        return 0
 
def regmul(A,B):
    res=np.matmul(A,B)
    return np.where(abs(res)<tol,0,res)
def main():
    mat=np.array(((1,0,0),(0,1,0),(0,1,0)),dtype=np.float)

    L,U,perms=LUdecomp(mat)
    mat=np.array([[1,0,0],[0,0,1],[0,1,0]],dtype=np.float)
    #mat=np.array([[1,1,0,2,3],[1,0,0,2,5]],dtype=np.float)

    mat=np.array([[2,0,1],[0,1,2]],dtype=np.float)
    mat=np.array([[1,1,2],[1,0,3]],dtype=np.float)
    b=np.transpose(np.array([[2,4],[6,0]],dtype=np.float))
    #vals,vecs=eigenproblem(mat+np.transpose(mat))
    Q,R,r=QRdecomp(np.transpose(mat))
    print(Q,"--\n",R,"--\n",r)
    print("\n___\n")
    print(np.matmul(Q,R))
    print(QRsolve(mat,b))

    """
    vals,vecs=eigenproblem(mat)
    print("eigenvalues: ",vals)
    print("eigenvectors:",vecs)

    Q,R,r=QRdecomp(mat)

    print(Q,"--\n",R,"--\n",r)
    print("\n___\n")
    print(regmul(Q,np.transpose(Q)))
    print(regmul(Q,R))

    print(QRsolve2(mat,b))

    Q,R,r=QRdecomp2(mat)
    print(mat)
    print("_________")
    print(Q)
    print("_________")
    print(R)
    print("_________")
    print(np.matmul(Q,R))
    input()
    print("rank:",r)

    #print(np.matmul(Q,np.transpose(Q)))
    #Q,R=QRdecomp(mat)

    imat=LUPinverse(mat)
    print(LUPinverse(mat))
    print(np.matmul(mat,imat))

    Q,R=QRdecomp(mat)
    print("mat:",mat)
    print("Q:",Q,R)
    print(np.matmul(Q,R))
    """
if __name__=="__main__":
    main()