import numpy as np
from matplotlib import pyplot as plt
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath

def M(a,b,z):
	N = 10
	M = 0
	for i in range(N):
		M += nFunc(a,i)/nFunc(b,i)*z**i/misc.factorial(i)
	return M
	
def nFunc(a,n):
	ans = 1
	for i in range(n):
		ans = ans*(a+i)
	return ans

def y1(a,xi):
	return np.exp(-0.25*xi**2)*M(0.5*a+0.25,0.5,0.5*xi**2)

def y2(a,xi):
	return np.exp(-0.25*xi**2)*xi*M(0.5*a+0.75,1.5,0.5*xi**2)
	
def Y1(a,xi):
	return y1(a,xi)+np.sqrt(-1j*1j*a)*y2(a,xi)
	
def Y2(a,xi):
	return y1(a,xi)-np.sqrt(-1j*1j*a)*y2(a,xi)	
	
def dY1(x,a,c,X):
	return Y1(a,c*(x+X))

def dY2(x,a,c,X):
	return Y2(a,c*(x+X))
	
def fun(E_,k,nu,delta,L,Z,phi):
	if isinstance(E_,float):
		E = E_
	else:
		E = E_[0]
    
	if abs(E)>abs(delta):
		return 1
	
	if np.isnan(E):
		print('E is nan')
		return 1
		
	a_e = -(E+nu/2)
	a_h = -(-E+nu/2)	
		
	xiL_e = np.sqrt(2/nu)*(0+k*nu) #or -k*nu*sgn?
	xiL_h = np.sqrt(2/nu)*(0-k*nu) #or +k*nu*sgn?
	xiR_e = np.sqrt(2/nu)*(L+k*nu) #or -k*nu*sgn?
	xiR_h = np.sqrt(2/nu)*(L-k*nu) #or +k*nu*sgn?
	
	pre_e1 = np.exp(0.5*k**2*nu-np.sqrt(-1j*1j*2*a_e*nu)*k)
	pre_e2 = np.exp(0.5*k**2*nu+np.sqrt(-1j*1j*2*a_e*nu)*k)
	pre_h1 = np.exp(0.5*k**2*nu+np.sqrt(-1j*1j*2*a_e*nu)*k)
	pre_h2 = np.exp(0.5*k**2*nu-np.sqrt(-1j*1j*2*a_e*nu)*k)
	
	D1_eL = pre_e1*Y1(a_e,xiL_e)
	D2_eL = pre_e2*Y2(a_e,xiL_e)
	D1_hL = pre_h1*Y1(a_h,xiL_h)
	D2_hL = pre_h2*Y2(a_h,xiL_h)
	
	D1_eR = pre_e1*Y1(a_e,xiR_e)
	D2_eR = pre_e2*Y2(a_e,xiR_e)
	D1_hR = pre_h1*Y1(a_h,xiR_h)
	D2_hR = pre_h2*Y2(a_h,xiR_h)
	
	#####################################
	kp=np.sqrt(1+2*E/nu)
	km=np.sqrt(1-2*E/nu)
	kx_e = kp+1j*k
	kx_h = km-1j*k
	#print('Should be 1: ',D1_eL)
	#print('Should be 1: ',D2_eL)
	#print('Should be',np.exp(1j*kx_e*L),": ",D1_eR)
	#print('Should be',np.exp(-1j*kx_e*L),": ",D2_eR)
	######################################
	
	dD1_eL = pre_e1*misc.derivative(dY1,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.001)
	dD2_eL = pre_e2*misc.derivative(dY2,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.001)
	dD1_hL = pre_h1*misc.derivative(dY1,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.001)
	dD2_hL = pre_h2*misc.derivative(dY2,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.001)
	
	dD1_eR = pre_e1*misc.derivative(dY1,L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.001)
	dD2_eR = pre_e2*misc.derivative(dY2,L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.001)
	dD1_hR = pre_h1*misc.derivative(dY1,L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.001)
	dD2_hR = pre_h2*misc.derivative(dY2,L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.001)
	
	ZL = Z
	ZR = Z
	phiL = 0
	phiR = phi
	
	qe = np.sqrt(1-k*k+1j*2/nu*np.sqrt(delta*delta-E*E))
	qh = np.sqrt(1-k*k-1j*2/nu*np.sqrt(delta*delta-E*E))

	eta = np.arccos(E/delta)

	gamma_eL = np.exp(1j*(-eta-phiL))
	gamma_eR = np.exp(1j*(-eta-phiR))
	gamma_hL = np.exp(1j*(eta-phiL))
	gamma_hR = np.exp(1j*(eta-phiR))
	   
	matrix = np.array([[D1_eL, D2_eL, 0, 0, gamma_hL, gamma_eL, 0, 0],
					   [0, 0, D1_hL, D2_hL, 1, 1, 0, 0],
					   [dD1_eL, dD2_eL, 0, 0, (ZL+1j)*qh*gamma_hL, (Z-1j)*qe*gamma_eL, 0, 0],
					   [0, 0, dD1_hL, dD2_hL, (ZL+1j)*qh, (Z-1j)*qe, 0, 0],
					   [D1_eR, D2_eR, 0, 0, 0, 0, gamma_hR, gamma_eR],
					   [0, 0, D1_hR, D2_hR, 0, 0, 1, 1],
					   [dD1_eR, dD2_eR, 0, 0, 0, 0, (-ZR-1j)*qh*gamma_hR, (-ZR+1j)*qe*gamma_eR],
					   [0, 0, dD1_hR, dD2_hR, 0, 0, (-ZR-1j)*qh, (-ZR+1j)*qe]])
	#print(matrix,'\n')
	D = det(matrix)
	return np.absolute(D)#D.real#D.real+D.imag

def makePlotPhi(nu, delta, Z, k, L, N, n):
	buf = 0.1
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	E_array = np.zeros(phi_array.shape)
	e0 = np.linspace(-delta+buf*delta,delta-buf*delta,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		
		for i in range(N):
			"""
			rootResult = opt.root(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2*delta #this is a quick fix
			"""
			E_array[i] = opt.fsolve(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))[0]
		print(E_array)
		plt.plot(phi_array,E_array,'.b')
	plt.axis([phi_start,phi_end,-delta,delta])
	title = 'nu = '+str(nu)
	plt.title(title)
	plt.show()

hw = 0.01#0.155789473684
N = 500
n = 3
nu = 100000.#40./hw
delta = 20.#2./hw
Z = 0
phi = 1.
L = 106.7
k = 0.

makePlotPhi(nu,delta,Z,k,L,N,n)