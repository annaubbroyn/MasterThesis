import numpy as np
from scipy import integrate
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
from scipy import optimize as opt
import mpmath
#import math

def y1(a,xi):
	res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
	factor = np.exp(-0.25*xi**2)
	return (float(factor.real*res.real-factor.imag*res.imag)+1j*float(factor.real*res.imag+factor.imag*res.real))
	
def y2(a,xi):
	res = xi*mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	factor = np.exp(-0.25*xi**2)
	return (float(factor.real*res.real-factor.imag*res.imag)+1j*float(factor.real*res.imag+factor.imag*res.real))
	
def Y1(a,xi,method):
	if method == 'y1y2':
		return y1(a,xi)
	elif method == 'pcfuToy1y2':
		Ua0 = np.sqrt(np.pi)/(mpmath.power(2,0.5*a+0.25)*mpmath.gamma(3./4.+0.5*a))
		C1 = (1-np.sin(np.pi*a))/(2*Ua0)
		C2 = 1/(2*Ua0)
		U = mpmath.pcfu(a,xi)
		Um = mpmath.pcfu(a,-xi)
		V = np.sin(np.pi*a)*U+Um
		return float(C1*U+C2*V)
	elif method == 'pcfu':
		return float(mpmath.pcfu(a,xi))
			
def Y2(a,xi,method):
	if method == 'y1y2':
		return y2(a,xi)
	elif method == 'pcfuToy1y2':
		dUa0 = -np.sqrt(np.pi)/(mpmath.power(2,0.5*a-0.25)*mpmath.gamma(0.25+0.5*a))
		C1 = (1+np.sin(np.pi*a))/(2*dUa0)
		C2 = -1/(2*dUa0)
		U = mpmath.pcfu(a,xi)
		Um = mpmath.pcfu(a,-xi)
		V = np.sin(np.pi*a)*U+Um
		return float(C1*U+C2*V)
	elif method == 'pcfu':
		if np.isinf(sp.gamma(1./2.+a)):
			return np.inf
		return float((mpmath.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*mpmath.pcfu(a,xi)+mpmath.pcfu(a,-xi)))
	
def dY1(x,a,c,X,method):
	return Y1(a,c*(x+X),method)

def dY2(x,a,c,X,method):
	return Y2(a,c*(x+X),method)

def fun(E_,phi,B,ky,Ef,L,Z,method):

	if isinstance(E_,float):
		E = E_
	else:
		E = E_[0]
	
	#if abs(E)>1:
	#	return [2*(abs(E)-1),2*(abs(E)-1)]
	
	
	if E>1:
		E = 1.
	elif E<-1:
		E = -1.
	
	nu = 2*Ef/B
	beta = E/Ef
	alpha = L/nu
	
	a_e = -nu/2*(1+beta)
	a_h = -nu/2*(1-beta)	
		
	xiL_e = np.sqrt(2*nu)*(0+ky)
	xiL_h = np.sqrt(2*nu)*(0-ky)
	xiR_e = np.sqrt(2*nu)*(alpha+ky)
	xiR_h = np.sqrt(2*nu)*(alpha-ky)
	
	D1_eL = Y1(a_e,xiL_e,method)
	D2_eL = Y2(a_e,xiL_e,method)
	D1_hL = Y1(a_h,xiL_h,method)
	D2_hL = Y2(a_h,xiL_h,method)
	
	D1_eR = Y1(a_e,xiR_e,method)
	D2_eR = Y2(a_e,xiR_e,method)
	D1_hR = Y1(a_h,xiR_h,method)
	D2_hR = Y2(a_h,xiR_h,method)
	
	"""
	dD1_eL = (1/nu)*misc.derivative(dY1,0,args=(a_e,np.sqrt(2*nu),ky,method),dx=0.001)
	dD2_eL = (1/nu)*misc.derivative(dY2,0,args=(a_e,np.sqrt(2*nu),ky,method),dx=0.001)
	dD1_hL = (1/nu)*misc.derivative(dY1,0,args=(a_h,np.sqrt(2*nu),-ky,method),dx=0.001)
	dD2_hL = (1/nu)*misc.derivative(dY2,0,args=(a_h,np.sqrt(2*nu),-ky,method),dx=0.001)
	
	dD1_eR = (1/nu)*misc.derivative(dY1,alpha,args=(a_e,np.sqrt(2*nu),ky,method),dx=0.001)
	dD2_eR = (1/nu)*misc.derivative(dY2,alpha,args=(a_e,np.sqrt(2*nu),ky,method),dx=0.001)
	dD1_hR = (1/nu)*misc.derivative(dY1,alpha,args=(a_h,np.sqrt(2*nu),-ky,method),dx=0.001)
	dD2_hR = (1/nu)*misc.derivative(dY2,alpha,args=(a_h,np.sqrt(2*nu),-ky,method),dx=0.001)
	"""
	"""
	dD1_eL = misc.derivative(dY1,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*ky,method),dx=0.001)
	dD2_eL = misc.derivative(dY2,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*ky,method),dx=0.001)
	dD1_hL = misc.derivative(dY1,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*ky,method),dx=0.001)
	dD2_hL = misc.derivative(dY2,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*ky,method),dx=0.001)
	
	dD1_eR = misc.derivative(dY1,L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*ky,method),dx=0.001)
	dD2_eR = misc.derivative(dY2,L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*ky,method),dx=0.001)
	dD1_hR = misc.derivative(dY1,L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*ky,method),dx=0.001)
	dD2_hR = misc.derivative(dY2,L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*ky,method),dx=0.001)
	"""
	dD1_eL = np.sqrt(2/nu)*misc.derivative(dY1,xiL_e,args=(a_e,1,0,method),dx=0.001)
	dD2_eL = np.sqrt(2/nu)*misc.derivative(dY2,xiL_e,args=(a_e,1,0,method),dx=0.001)
	dD1_hL = np.sqrt(2/nu)*misc.derivative(dY1,xiL_h,args=(a_h,1,0,method),dx=0.001)
	dD2_hL = np.sqrt(2/nu)*misc.derivative(dY2,xiL_h,args=(a_h,1,0,method),dx=0.001)
	
	dD1_eR = np.sqrt(2/nu)*misc.derivative(dY1,xiR_e,args=(a_e,1,0,method),dx=0.001)
	dD2_eR = np.sqrt(2/nu)*misc.derivative(dY2,xiR_e,args=(a_e,1,0,method),dx=0.001)
	dD1_hR = np.sqrt(2/nu)*misc.derivative(dY1,xiR_h,args=(a_h,1,0,method),dx=0.001)
	dD2_hR = np.sqrt(2/nu)*misc.derivative(dY2,xiR_h,args=(a_h,1,0,method),dx=0.001)
	
	ZL = Z
	ZR = Z
	phiL = 0
	phiR = phi
	
	qe = np.sqrt(1-ky*ky+1j*np.sqrt(1-E**2)/Ef)
	qh = np.sqrt(1-ky*ky-1j*np.sqrt(1-E**2)/Ef)
	
	eta = np.arccos(E)

	gamma_eL = np.exp(1j*(-eta-phiL))
	gamma_eR = np.exp(1j*(-eta-phiR))
	gamma_hL = np.exp(1j*(eta-phiL))
	gamma_hR = np.exp(1j*(eta-phiR))
	
	   
	matrix = np.array([[D1_eL, D2_eL, 0, 0, gamma_hL, gamma_eL, 0, 0],
					   [0, 0, D1_hL, D2_hL, 1, 1, 0, 0],
					   [dD1_eL, dD2_eL, 0, 0, (ZL+1j)*qh*gamma_hL, (ZL-1j)*qe*gamma_eL, 0, 0],
					   [0, 0, dD1_hL, dD2_hL, (ZL+1j)*qh, (ZL-1j)*qe, 0, 0],
					   [D1_eR, D2_eR, 0, 0, 0, 0, gamma_hR, gamma_eR],
					   [0, 0, D1_hR, D2_hR, 0, 0, 1, 1],
					   [dD1_eR, dD2_eR, 0, 0, 0, 0, (-ZR-1j)*qh*gamma_hR, (-ZR+1j)*qe*gamma_eR],
					   [0, 0, dD1_hR, dD2_hR, 0, 0, (-ZR-1j)*qh, (-ZR+1j)*qe]])
	
	D = det(matrix)
	
	return [D.real,D.imag]

def freeEnergy(phi,ky,B,Ef,L,Z,kBT,method):
	n=4
	de0 = 0.01
	e0 = np.linspace(-1+de0,1-de0,n)
	temp_E_array = np.zeros(n)
	E_array = np.zeros(2)
	success = np.zeros(n,dtype=bool)
	maxDiff = 10**(-6)
	num_success = 0
	#for j in range(n):
	for j in range(n):
		rootResult = opt.root(fun,e0[j],args=(phi,B,ky,Ef,L,Z,method))
		temp_E_array[j] = rootResult.x[0]
		success[j] = rootResult.success
		if rootResult.success:
			num_success += 1
	if num_success == n:
		E_array[0] = temp_E_array[0]
		E_array[1] = temp_E_array[-1]
	elif num_success == 0:
		E_array[0] = -1.
		E_array[1] = 1.
	else:
		nextIndex = 0
		for j in range(n):
			if success[j]:
				E_array[0] = temp_E_array[j]
				index = j
				break
		for j in range(index,n):
			if (success[j] and np.absolute(E_array[0]-temp_E_array[j])>10**(-6)):
				E_array[1] = temp_E_array[j]
				break
			if j == n-1:
				if temp_E_array[0]>0.9999:
					E_array[1] = 1.
				elif temp_E_array[0]<-0.9999:
					E_array[1] = -1.
				else:
					E_array[1] = E_array[0]
	E_array[0] = min(E_array[0],1.)
	E_array[0] = max(E_array[0],-1.)
	E_array[1] = min(E_array[1],1.)
	E_array[1] = max(E_array[1],-1.)
			
	result = 0
	for i in range(2):
		result += -np.log(2*np.cosh(0.5*E_array[i]/kBT))
	return result					
	
	
def dFreeEnergy(ky,phi,B,Ef,L,Z,kBT,method):
	return misc.derivative(freeEnergy,phi,args=(ky,B,Ef,L,Z,kBT,method),dx=0.001)
	
def totalCurrent(phi,B,Ef,L,Z,kBT,method):
	kyMin = -0.01
	kyMax = 0.01
	return integrate.quad(dFreeEnergy,kyMin,kyMax,args=(phi,B,Ef,L,Z,kBT,method))[0]
