import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath

"""
def tempU(a,x):
    return float(mpmath.pcfd(-a-1./2.,x))
	
def tempV(a,x):
    if np.isinf(sp.gamma(1./2.+a)):
        return np.inf
    return (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*tempU(a,x)+tempU(a,-x))

def U(a,xi):
#Actually y1
	Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	C1 = (1-np.sin(np.pi*a))/(2*Ua0)
	C2 = np.pi/(sp.gamma(0.5+a)*2*Ua0)
	return C1*tempU(a,xi)+C2*tempV(a,xi)
	
def V(a,xi):
#Actually y2
	Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	C1 = (1+np.sin(np.pi*a))/(2*dUa0)
	C2 = -np.pi/(sp.gamma(0.5+a)*2*dUa0)
	return C1*U(a,xi)+C2*tempV(a,xi)
"""
##################################

def U(a,xi):
#Actually my y1
	return np.exp(-0.25*xi**2)*myM(0.5*a+0.25,0.5,0.5*xi**2)

def V(a,xi):
#Actually my y2
	return np.exp(-0.25*xi**2)*xi*myM(0.5*a+0.75,1.5,0.5*xi**2)


##################################

def myM(a,b,z):
	N = 50
	M = 0
	for i in range(N):
		M += nFunc(a,i)/nFunc(b,i)*z**i/misc.factorial(i)
	return M
	
def nFunc(a,n):
	ans = 1
	for i in range(n):
		ans = ans*(a+i)
	return ans
"""
##################################

def U(a,z):
	U0a = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	dU0a = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	return np.exp(-0.25*z**2)*(U0a*myM(0.5*a+0.25,0.5,0.5*z**2)+dU0a*z*myM(0.5*a+3./4.,3./2.,0.5*z**2))


"""
###############################
"""
def U(a,x):
   return float(mpmath.pcfd(-a-1./2.,x))

def V(a,x):
    if np.isinf(sp.gamma(1./2.+a)):
        return np.inf
    return (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U(a,x)+U(a,-x))
"""
def dUFunction(x,a,c,X):
    return U(a,c*(x+X))

def dVFunction(x,a,c,X):
    return V(a,c*(x+X))
	
def fun(E_,k,nu,delta,L,Z,phi):
	if isinstance(E_,float):
		E = E_
	else:
		E = E_[0]
    
	if abs(E)>abs(delta):
		return 1
		
	a_e = -(E+nu/2)
	a_h = -(-E+nu/2)	
		
	xiL_e = np.sqrt(2/nu)*(0+k*nu) #or -k*nu*sgn?
	xiL_h = np.sqrt(2/nu)*(0-k*nu) #or +k*nu*sgn?
	xiR_e = np.sqrt(2/nu)*(L+k*nu) #or -k*nu*sgn?
	xiR_h = np.sqrt(2/nu)*(L-k*nu) #or +k*nu*sgn?
	
	D1_eL = U(a_e,xiL_e)
	D2_eL = V(a_e,xiL_e)
	D1_hL = U(a_h,xiL_h)
	D2_hL = V(a_h,xiL_h)
	
	D1_eR = U(a_e,xiR_e)
	D2_eR = V(a_e,xiR_e)
	D1_hR = U(a_h,xiR_h)
	D2_hR = V(a_h,xiR_h)
	
	dD1_eL = misc.derivative(dUFunction,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.0001)
	dD2_eL = misc.derivative(dVFunction,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.0001)
	dD1_hL = misc.derivative(dUFunction,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.0001)
	dD2_hL = misc.derivative(dVFunction,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.0001)
	
	ZL = Z
	ZR = Z

	eta = np.arccos(E/delta)

	gamma_eL = np.exp(-1j*eta)
	gamma_hL = np.exp(1j*eta)

	qe = np.sqrt(1-k*k+1j*2/nu*np.sqrt(delta*delta-E*E))
	qh = np.sqrt(1-k*k-1j*2/nu*np.sqrt(delta*delta-E*E))

	matrix = np.array([[D1_eL, D2_eL, 0, 0, gamma_hL, gamma_eL],
					   [0, 0, D1_hL, D2_hL, 1, 1],
					   [dD1_eL, dD2_eL, 0, 0, (ZL+1j)*qh*gamma_hL, (Z-1j)*qe*gamma_eL],
					   [0, 0, dD1_hL, dD2_hL, (ZL+1j)*qh,(Z-1j)*qe],
					   [D1_eR, D2_eR, 0, 0, 0, 0],
					   [0, 0, D1_hR, D2_hR, 0, 0]])

	D = det(matrix)
	return np.absolute(D)
		
def makePlotk(nu, delta, Z, phi, L, N, n):
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
			rootResult = opt.root(fun,e0[j],args=(k_array[i],nu,delta,L,Z,phi))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2*delta #this is a quick fix			
			#E_array[i] = opt.fsolve(fun,e0[j],args=(k_array[i],nu,delta,L,Z,phi))[0]
		plt.plot(k_array,E_array,'.b')
	plt.axis([k_start,k_end,0.05,delta-0.05])	
	plt.show()

	
	
N = 500
n = 4
nu = 40.#500.#40.#160.
delta = 2.#8.
Z = 0
phi = 1.
L = 106.7
k = 0.2

print('Solve_SNV with lib y1 and y2')
makePlotk(nu,delta,Z,phi,L,N,n)







   


