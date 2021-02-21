## transposée de A: second membre du pb adjoint
def f_dual(i,u_i,ystate,p_i):
    a11=-u_i*ystate[i,1]/N
    a12=-u_i*ystate[i,0]/N
    a21=u_i*ystate[i,1]/N
    a22=u_i*ystate[i,0]/N-beta
    A=np.matrix([[a11, a12, 0],[a21 ,a22, 0],[ 0 ,beta, 0]])
    return np.transpose(np.transpose(A)*np.transpose(np.matrix(p_i))+np.matrix([[0],[ystate[i,1]],[0]]))

def RK4_adj(u,ystate,NT):
    # solved in the time interval [a,b]
    # initial conditions = vector alpha
    # number of subdivisions in the interval [a,b] = N
    h = T/NT
    pT = np.zeros((NT,3))
    for i in list(reversed(range(1,NT))):
        k1 = h*f_dual(i,u[i],ystate,pT[i,:])
        k2 = h*f_dual(i,0.5*(u[i]+u[i-1]),ystate,pT[i,:]+0.5*k1)
        k3 = h*f_dual(i,0.5*(u[i]+u[i-1]),ystate,pT[i,:]+0.5*k2)
        k4 = h*f_dual(i,u[i-1],ystate,pT[i,:]+k3)
        pT[i-1,:] = pT[i,:] + (k1 + 2*k2 + 2*k3 + k4)/6

    return pT