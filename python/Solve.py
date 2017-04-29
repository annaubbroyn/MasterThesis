import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath

def U(a,x):
    denominator = 1#2**(-1./4.-0.5*a)*sp.gamma(1/4.-0.5*a)
    return float(mpmath.pcfd(-a-1./2.,x))

def dUFunction(x,a,c,X):
    return U(a,c*(x+X))

def V(a,x):
    if np.isinf(sp.gamma(1./2.+a)):
        return np.inf
    return (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U(a,x)+U(a,-x))

def dVFunction(x,a,c,X):
    return V(a,c*(x+X))
	
def fun(E_,k,nu,delta,W,Z,sgn):
	if isinstance(E_,float):
		E = E_
	else:
		E = E_[0]
    
	if abs(E)>abs(delta):
		return 1
	
	a_e = -(E+nu/2)
	a_h = -(-E+nu/2)
	eta = 2/nu*np.sqrt(delta*delta-E*E)
	gamma_e = delta/(E+1j*np.sqrt(delta*delta - E*E))
	gamma_h = delta/(E-1j*np.sqrt(delta*delta - E*E))
	q_e = np.sqrt(1-k*k-1j*eta)	# or -1j ?
	q_h = np.sqrt(1-k*k+1j*eta)	# or +1j ?
    
	xiW_e = np.sqrt(2/nu)*(W+k*nu*sgn) #or -k*nu*sgn?
	xiW_h = np.sqrt(2/nu)*(W-k*nu*sgn) #or +k*nu*sgn?
	xi0_e = np.sqrt(2/nu)*(0+k*nu*sgn) #or -k*nu*sgn?
	xi0_h = np.sqrt(2/nu)*(0-k*nu*sgn) #or +k*nu*sgn?

	U_e = U(a_e,xi0_e)
	U_h = U(a_h,xi0_h)
	dU_e = misc.derivative(dUFunction,0,args=(a_e,np.sqrt(2/nu),nu*k*sgn))
	dU_h = misc.derivative(dUFunction,0,args=(a_h,np.sqrt(2/nu),-nu*k*sgn))
	V_e = V(a_e,xi0_e)
	V_h = V(a_h,xi0_h)
	dV_e = misc.derivative(dVFunction,0,args=(a_e,np.sqrt(2/nu),nu*k*sgn))
	dV_h = misc.derivative(dVFunction,0,args=(a_h,np.sqrt(2/nu),-nu*k*sgn))
	

	Vfactor_e = 0    
	if np.isinf(V(a_e,xiW_e)):
		Vfactor_e = 0
		V_e=0
		dV_e = 0
	elif abs(U(a_e,xiW_e)) > 0 and abs(V(a_e,xiW_e)) > 0:
		Vfactor_e = U(a_e,xiW_e)/V(a_e,xiW_e)
		if Vfactor_e == 0.:
			V_e = 0
			dV_e = 0
	else:
		print('Error: V_W is zero but not U_W for electron')
    
	Vfactor_h = 0
	if np.isinf(V(a_h,xiW_h)):
		Vfactor_h = 0
		V_h =  0
		dV_h = 0
	elif abs(U(a_h,xiW_h)) > 0 and abs(V(a_h,xiW_h)) > 0:
		Vfactor_h = U(a_h,xiW_h)/V(a_h,xiW_h)
		if Vfactor_h == 0.:
			V_h = 0
			dV_h = 0
	else:
		print('Error: V_W is zero but not U_W for hole')
	
	chi_e = U_e-Vfactor_e*V_e
	chi_h = U_h-Vfactor_h*V_h
	dchi_e = dU_e-Vfactor_e*dV_e
	dchi_h = dU_h-Vfactor_h*dV_h
	
	mat_e = np.array([[chi_e,1],[dchi_e,Z+1j*q_e]])
	mat_h = np.array([[chi_h, 1],[dchi_h,Z-1j*q_h]])
	
	De = det(mat_e)
	Dh = det(mat_h)

	return (gamma_e*De*Dh).imag
		
def makePlot(nu, delta, Z, sgn, W, N, n):
	buf = 0.25
	k_start = -1.25
	k_end = 1.25
	k_array = np.linspace(k_start, k_end, N)
	E_array = np.zeros(k_array.shape)
	e0 = np.linspace(0+buf,delta-buf,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):		
			E_array[i] = opt.fsolve(fun,e0[j],args=(k_array[i],nu,delta,W,Z,sgn))[0]
			plt.plot(k_array,E_array,'.b')
	plt.axis([k_start,k_end,0.05,delta-0.05])	
	plt.show()

N = 500
n = 3
nu = 40.#160.
delta = 2.#8.
Z = 0
sgn = 1
W = 106.7

makePlot(nu,delta,Z,sgn,W,N,n)

"""

N = 100
k = -0.1
delta = 2.
nu = 40
Z = 0
W = 0
sgn = 1

E = np.linspace(0.,delta,N)
res = np.zeros(E.shape)

for i in range(N):
    res[i] = fun(E[i], k, nu, delta, W, Z, sgn)

plt.plot(E,res)
plt.show()

"""






   


