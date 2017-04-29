import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath
	
def fun(E_,k,nu,delta,L,Z,phi):
	if isinstance(E_,float):
		E = E_
	else:
		E = E_[0]
    
	if abs(E)>abs(delta):
		return 1
	
	ZL = Z
	ZR = Z
	phiL = 0
	phiR = phi
	
	kx = np.sqrt(1+1j*1j*k*k)
	kx_e = kx#np.sqrt(1+1j*1j*k*k+2*E/nu)
	kx_h = kx#np.sqrt(1+1j*1j*k*k-2*E/nu)
	
	
	ekR1_e = np.exp(1j*kx_e*L)
	ekR2_e = np.exp(-1j*kx_e*L)
	ekR1_h = np.exp(1j*kx_h*L)
	ekR2_h = np.exp(-1j*kx_h*L)
	
	qe = np.sqrt(1-k*k+1j*2/nu*np.sqrt(delta*delta-E*E))
	qh = np.sqrt(1-k*k-1j*2/nu*np.sqrt(delta*delta-E*E))
	
	eta = np.arccos(E/delta)
	
	gamma_eL = np.exp(1j*(-eta-phiL))
	gamma_eR = np.exp(1j*(-eta-phiR))
	gamma_hL = np.exp(1j*(eta-phiL))
	gamma_hR = np.exp(1j*(eta-phiR))
	
	matrix = np.array([[1, 1, 0, 0, gamma_hL, gamma_eL, 0, 0],
					   [0, 0, 1, 1, 1, 1, 0, 0],
					   [1j*kx_e, -1j*kx_e, 0, 0, (ZL+1j)*qh*gamma_hL, (Z-1j)*qe*gamma_eL, 0, 0],
					   [0, 0, 1j*kx_h, -1j*kx_h, (ZL+1j)*qh, (Z-1j)*qe, 0, 0],
					   [ekR1_e, ekR2_e, 0, 0, 0, 0, gamma_hR, gamma_eR],
					   [0, 0, ekR1_h, ekR2_h, 0, 0, 1, 1],
					   [1j*kx_e*ekR1_e, -1j*kx_e*ekR2_e, 0, 0, 0, 0, (-ZR-1j)*qh*gamma_hR, (-ZR+1j)*qe*gamma_eR],
					   [0, 0, 1j*kx_h*ekR1_h, -1j*kx_h*ekR2_h, 0, 0, (-ZR-1j)*qh, (-ZR+1j)*qe]])
	"""
	matrix = np.array([[1, 1, 0, 0, gamma_hL, gamma_eL, 0, 0],
					   [0, 0, 1, 1, 1, 1, 0, 0],
					   [1j*kx, -1j*kx, 0, 0, (ZL+1j)*qh*gamma_hL, (Z-1j)*qe*gamma_eL, 0, 0],
					   [0, 0, 1j*kx, -1j*kx, (ZL+1j)*qh, (Z-1j)*qe, 0, 0],
					   [ekR1, ekR2, 0, 0, 0, 0, gamma_hR, gamma_eR],
					   [0, 0, ekR1, ekR2, 0, 0, 1, 1],
					   [1j*kx*ekR1, -1j*kx*ekR2, 0, 0, 0, 0, (-ZR-1j)*qh*gamma_hR, (-ZR+1j)*qe*gamma_eR],
					   [0, 0, 1j*kx*ekR1, -1j*kx*ekR2, 0, 0, (-ZR-1j)*qh, (-ZR+1j)*qe]])
	"""
	D = det(matrix)
	return np.absolute(D)

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
			rootResult = opt.root(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2*delta #this is a quick fix
		plt.plot(phi_array,E_array,'.b')
	plt.axis([phi_start,phi_end,-delta-0.01*delta,delta+0.01*delta])	
	plt.show()


def testFunction(nu, delta, Z, k, L, N, n):
	phi_start = -np.pi+0.09
	phi_end = -np.pi+0.1
	phi = np.linspace(phi_start,phi_end,n)
	E_array = np.linspace(-delta,delta,N)
	Func_array = np.zeros(E_array.shape)
	
	plt.figure()
	for j in range(n):
		for i in range(N):		
			Func_array[i] = fun(E_array[i],k,nu,delta,L,Z,phi[j])
		plt.plot(E_array,Func_array,'.',label='phi='+str(phi[j]))
		plt.legend()
	plt.show()
	
hw = 0.01#0.155789473684
N = 500
n = 3
nu = 40./hw
delta = 2.#./hw
Z = 0
phi = 1.
L = 106.7
k = 0.2

makePlotPhi(nu,delta,Z,k,L,N,n)
#testFunction(nu,delta,Z,k,L,N,n)





   


