import numpy as np


def dirichlet(j,A,B,r,dx,dt,Nx,g1,g2):
    
    A[0,0] = 1-2*r
    A[0,1] = r
    A[-1,-1] = 1-2*r 
    A[-1,-2] = r
    
    for i in range(Nx+1):
        
        #dirichlet
        if(i==0):
            B[i]= r*g1(j*dt)
        
        #dirichlet
        elif(i==Nx):
            B[i] = r*g2(j*dt)
    
        else:
            B[i] = 0
        
        return None
    
def neumann(j,A,B,r,dx,dt,Nx,g1,g2):
    
    A[0,0] = 1-r
    A[0,1] = r
    A[-1,-1] = 1-r 
    A[-1,-2] = r
    
    for i in range(Nx+1):
        
        #Neumann
        if(i==0):
            B[i]= r*dx*g1(j*dt)
        
        #Neumann
        elif(i==Nx):
            B[i] = r*dx*g2(j*dt)
    
        else:
            B[i] = 0
        
        return None
 
    
def dirineumann(j,A,B,r,dx,dt,Nx,g1,g2):
    
    A[0,0] = 1-2*r
    A[0,1] = r
    A[-1,-1] = 1-r 
    A[-1,-2] = r
    
    for i in range(Nx+1):
        
        #Dirichlet
        if(i==0):
            B[i]= r*g1(j*dt)
        
        #Neumann
        elif(i==Nx):
            B[i] = r*dx*g2(j*dt)
    
        else:
            B[i] = 0
        
        return None
    

def FTCS(x0,xn,dx,t0,tn,dt,fini,g1,g2,bc):
    
    beta = 1
    Nx = int((xn-x0)/dx)
    x = np.linspace(x0,xn,Nx+1)
                   
    Nt = int((tn-t0)/dt)
    t = np.linspace(t0,tn,Nt+1)
    
    r = beta * dt/(dx)**2
    
    if(r>0.5):
        print(" r is greater than 1/2")
    
    U = np.zeros((Nt+1,Nx+1))
    
    #initial condition
    for i in range(Nx+1):
        U[0,i] = fini(x[i])
    
    #boundary condition:
    B = np.zeros(Nx+1)
    A = np.zeros((Nx+1,Nx+1))
    
    for i in range(1,Nx):
        A[i,i-1] = r
        A[i,i] = 1-2*r
        A[i,i+1] = r
    
    for j in range(Nt):
       
        if (bc=="DD"):
            dirichlet(j,A,B,r,dx,dt,Nx,g1,g2)
        elif(bc=="DN"):
            dirineumann(j,A,B,r,dx,dt,Nx,g1,g2)
        elif(bc=="NN"):
            neumann(j,A,B,r,dx,dt,Nx,g1,g2)
        else:
            print("wrong input")
            
        Uk = np.zeros(Nx+1)
        
        for i in range(Nx+1):
            Uk[i] = U[j,i]
        
        Uk1 = np.dot(A,Uk) + B
        
        for i in range(Nx+1):
            U[j+1,i] = Uk1[i]
            

    return t,x,U
            
    
    
        
    
    
    
    
    
    
    