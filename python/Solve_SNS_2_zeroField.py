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
    return float(mpmath.pcfd(-a-1./2.,x))

def dUFunction(x,a,c,X):
    return U(a,c*(x+X))

def V(a,x):
    if np.isinf(sp.gamma(1./2.+a)):
        return np.inf
    return (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U(a,x)+U(a,-x))

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
	
	kp=np.sqrt(1+2*E/nu)
	km=np.sqrt(1-2*E/nu)
	kx = np.sqrt(1+1j*1j*k*k)
	kx_e = kp+1j*k#np.sqrt(1+1j*1j*k*k+2*E/nu)#kx
	kx_h = km-1j*k#np.sqrt(1+1j*1j*k*k-2*E/nu)#kx
		
	xiL_e = np.sqrt(2/nu)*(0+k*nu) #or -k*nu*sgn?
	xiL_h = np.sqrt(2/nu)*(0-k*nu) #or +k*nu*sgn?
	xiR_e = np.sqrt(2/nu)*(L+k*nu) #or -k*nu*sgn?
	xiR_h = np.sqrt(2/nu)*(L-k*nu) #or +k*nu*sgn?
	
	D1_eL = 1
	D2_eL = 1
	D1_hL = 1
	D2_hL = 1
	
	D1_eR = np.exp(1j*kx_e*L)
	D2_eR = np.exp(-1j*kx_e*L)
	D1_hR = np.exp(-1j*kx_h*L)
	D2_hR = np.exp(1j*kx_h*L)
	
	dD1_eL = 1j*kx_e
	dD2_eL = -1j*kx_e
	dD1_hL = -1j*kx_h
	dD2_hL = 1j*kx_h
	
	dD1_eR = 1j*kx_e*np.exp(1j*kx_e*L)
	dD2_eR = -1j*kx_e*np.exp(-1j*kx_e*L)
	dD1_hR = -1j*kx_h*np.exp(-1j*kx_h*L)
	dD2_hR = 1j*kx_h*np.exp(1j*kx_h*L)
	
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
	D = det(matrix)
	return D.real+D.imag
		
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
			E_array[i] = opt.fsolve(fun,e0[j],args=(k_array[i],nu,delta,L,Z,phi))[0]
			plt.plot(k_array,E_array,'.b')
	plt.axis([k_start,k_end,0.05,delta-0.05])	
	plt.show()

def makePlotPhi(nu, delta, Z, k, L, N, n):
	buf = 0.1
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	E_array = np.zeros(phi_array.shape)
	e0 = np.linspace(-delta+buf,delta-buf,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		
		for i in range(N):
			rootResult = opt.root(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2*delta #this is a quick fix
		
		#	E_array[i] = opt.fsolve(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))[0]
		
		plt.plot(phi_array,E_array,'.b')
	plt.axis([phi_start,phi_end,-delta,delta])
	title = 'nu = '+str(nu)
	plt.title(title)
	plt.show()

def testFunction(nu, delta, Z, k, L, N, n):
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi = np.linspace(phi_start,phi_end,n)
	E_array = np.linspace(-delta,delta,N)
	Func_array = np.zeros(E_array.shape)
	plt.figure()
	for j in range(n):
		for i in range(N):		
			Func_array[i] = fun(E_array[i],k,nu,delta,L,Z,phi[j])
		plt.plot(E_array,Func_array)
	plt.show()
	
def plotAndSave(hw_start, hw_end, figCount, Z, k, L, N, n):
	hw_array = np.linspace(hw_start,hw_end,figCount)
	for hw in hw_array:
		print('hw:',hw)
		nu = 40./hw
		delta = 2./hw
		buf = 0.1
		phi_start = -3*np.pi
		phi_end = 3*np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		E_array = np.zeros(phi_array.shape)
		e0 = np.linspace(-delta+buf,delta-buf,n)
		fig = plt.figure()
		for j in range(n):
			print('begin:',j+1,'/',n)
			for i in range(N):		
				E_array[i] = opt.fsolve(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))[0]
			plt.plot(phi_array,E_array,'.b')
			print('finish:',j+1,'/',n)
		plt.axis([phi_start,phi_end,-delta,delta])
		title = 'nu = '+str(nu)
		plt.title(title)
		path = 'figures/overNight/'
		name = 'hw_'+str(hw)
		fig.savefig(path+name+'.png')
	

hw = 0.01#155789473684
N = 100
n = 10
nu = 40./hw
delta = 20.#/hw
Z = 0
phi = 1.
L = 106.7
k = 0.2

hw_start = 0.11
hw_end = 0.4
figCount = 20

print('hw:',hw)
#makePlotk(nu,delta,Z,phi,L,N,n)
makePlotPhi(nu,delta,Z,k,L,N,n)
#testFunction(nu, delta, Z, k, L, N, n)
#plotAndSave(hw_start, hw_end, figCount, Z, k, L, N, n)




   


