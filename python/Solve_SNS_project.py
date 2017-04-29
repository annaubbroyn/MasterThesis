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
	
	k = 1
	L = 1
	ZL = Z
	ZR = Z
	
	a_e = -(E+nu/2)
	a_h = -(-E+nu/2)
	eta = np.arccos(E/delta)

	gamma_eL = np.exp(1j*eta)
	gamma_eR = np.exp(1j*(eta-phi))
	gamma_hL = np.exp(1j*eta)
	gamma_hR = np.exp(1j*(eta+phi))

	ekRe = np.exp(1j*k*L)
	ekRh = np.exp(-1j*k*L)
	
	q_eL = np.sqrt(1-1j*2/nu*np.sqrt(delta*delta-E*E))	# or -1j ?
	q_hL = np.sqrt(1+1j*2/nu*np.sqrt(delta*delta-E*E))	# or +1j ?
	q_eR = -q_eL
	q_hR = -q_hL
	
	matrix = np.array([[1, 1, 0, 0, gamma_eL, 1, 0, 0],
					   [0, 0, 1, 1, 1, gamma_hL, 0, 0],
					   [1j*k, -1j*k, 0, 0, (ZL-1j)*k*gamma_eL, (Z+1j)*k, 0, 0],
					   [0, 0, 1j*k, -1j*k, (ZL-1j)*k, (Z+1j)*k*gamma_hL, 0, 0],
					   [ekRe, ekRh, 0, 0, 0, 0, gamma_eR, 1],
					   [0, 0, ekRe, ekRh, 0, 0, 1, gamma_hR],
					   [1j*k*ekRe, -1j*k*ekRh, 0, 0, 0, 0, (-ZR+1j)*k*gamma_eR, (-ZR-1j)*k],
					   [0, 0, 1j*k*ekRe, -1j*k*ekRh, 0, 0, (-ZR+1j)*k, (-ZR-1j)*k*gamma_hR]])
	"""
	matrix = np.array([[gamma_eL, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
					   [0, gamma_eL, 0, -1, 0, 0, 0, 0,0, 0, 1, 1, 0, 0, 0, 0],
					   [0, -1, 0, gamma_hL, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
					   [1, 0, gamma_hL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
					   [0, 0, 0, 0, gamma_eR, 0, 1, 0, ekRe, ekRh, 0, 0, 0, 0, 0, 0],
					   [0, 0, 0, 0, 0, gamma_eR, 0, -1, 0, 0, ekRe, ekRh, 0, 0, 0, 0],
					   [0, 0, 0, 0, 0, -1, 0, gamma_hR, 0, 0, 0, 0, ekRe, ekRh, 0, 0],
					   [0, 0, 0, 0, 1, 0, gamma_hR, 0, 0, 0, 0, 0, 0, 0, ekRe, ekRh],
					   [(-1j+ZL)*k*gamma_eL, 0, (1j+ZL)*k, 0, 0, 0, 0, 0, 1j*k, -1j*k, 0, 0, 0, 0, 0, 0],
					   [0, (-1j+ZL)*k*gamma_eL, 0, -(1j+ZL)*k, 0, 0, 0, 0,0, 0, 1j*k, -1j*k, 0, 0, 0, 0],
					   [0, -(-1j+ZL)*k, 0, (1j+ZL)*k*gamma_hL, 0, 0, 0, 0, 0, 0, 0, 0, 1j*k, -1j*k, 0, 0],
					   [(-1j+ZL)*k, 0, (1j+ZL)*k*gamma_hL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1j*k, -1j*k],
					   [0, 0, 0, 0, (1j-ZR)*k*gamma_eR, 0, (-1j-ZR)*k, 0, 1j*k*ekRe, -1j*k*ekRh, 0, 0, 0, 0, 0, 0],
					   [0, 0, 0, 0, 0, (1j-ZR)*k*gamma_eR, 0, -(-1j-ZR)*k, 0, 0, 1j*k*ekRe, -1j*k*ekRh, 0, 0, 0, 0],
					   [0, 0, 0, 0, 0, -(1j-ZR)*k, 0, (-1j-ZR)*k*gamma_hR, 0, 0, 0, 0, 1j*k*ekRe, -1j*k*ekRh, 0, 0],
					   [0, 0, 0, 0, (1j-ZR)*k, 0, (-1j-ZR)*k*gamma_hR, 0, 0, 0, 0, 0, 0, 0, 1j*k*ekRe, -1j*k*ekRh]])
	"""
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
			E_array[i] = opt.fsolve(fun,e0[j],args=(k_array[i],nu,delta,L,Z,phi))[0]
			plt.plot(k_array,E_array,'.b')
	plt.axis([k_start,k_end,0.05,delta-0.05])	
	plt.show()

def makePlotPhi(nu, delta, Z, k, L, N, n):
	buf = 0.25
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	E_array = np.zeros(phi_array.shape)
	e0 = np.linspace(-delta+buf,delta-buf,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):		
			E_array[i] = opt.fsolve(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))[0]
		plt.plot(phi_array,E_array,'.b')
	plt.axis([phi_start,phi_end,-delta,delta])	
	plt.show()

	
	
N = 300
n = 3
nu = 40.#160.
delta = 2.#8.
Z = 0.5
phi = 1.
L = 106.7
k = 1.3

#makePlotk(nu,delta,Z,phi,L,N,n)
makePlotPhi(nu,delta,Z,k,L,N,n)

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






   


