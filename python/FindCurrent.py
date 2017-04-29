import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy import integrate
from scipy.linalg import det
import mpmath

def y1(a,xi):
	res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
	return np.exp(-0.25*xi**2)*(float(res.real)+1j*float(res.imag))
	
def y2(a,xi):
	res = xi*mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	return np.exp(-0.25*xi**2)*(float(res.real)+1j*float(res.imag))
	
def Y1(a,xi):
	#Low field:
	#return y1(a,xi)+np.sqrt(-1j*1j*a)*y2(a,xi)
	
	#Only y1 (fungerer for nu>100):
	return y1(a,xi)
	
	#y1 from pcfd
	#Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	#dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	#C1 = (1-np.sin(np.pi*a))/(2*Ua0)
	#C2 = np.pi/(sp.gamma(0.5+a)*2*Ua0)
	#U = float(mpmath.pcfd(-a-1./2.,xi))
	#Um = float(mpmath.pcfd(-a-1./2.,xi))
	#V = (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U+Um)
	#return C1*U+C2*V
	
	#High field:
	#Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	#dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	#return Ua0*y1(a,xi)+dUa0*y2(a,xi)
	
	#Pcfd
	#return float(mpmath.pcfd(-a-1./2.,xi))
	
def Y2(a,xi):
	#Low field:
	#return y1(a,xi)-np.sqrt(-1j*1j*a)*y2(a,xi)	
	
	#Only y2:
	return y2(a,xi)
	
	#y2 from pcfd
	#Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	#dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	#C1 = (1+np.sin(np.pi*a))/(2*dUa0)
	#C2 = -np.pi/(sp.gamma(0.5+a)*2*dUa0)
	#U = float(mpmath.pcfd(-a-1./2.,xi))
	#Um = float(mpmath.pcfd(-a-1./2.,xi))
	#V = (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U+Um)
	#return C1*U+C2*V
	
	#High field:
	#Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	#dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	#return sp.gamma(0.5+a)/np.pi*((np.sin(np.pi*a)+1)*Ua0*y1(a,xi)+(np.sin(np.pi*a)-1)*dUa0*y2(a,xi))
	
	#Pcdf
	#if np.isinf(sp.gamma(1./2.+a)):
	#	return np.inf
	#return (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*Y1(a,xi)+Y1(a,-xi))

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
	
	if np.isnan(E) or np.isinf(E):
		print('E is ',E)
		return 1
	
	a_e = -(E+nu/2)
	a_h = -(-E+nu/2)	
	
	kp=np.sqrt(1+2*E/nu)
	km=np.sqrt(1-2*E/nu)
	kx = np.sqrt(1+1j*1j*k*k)
	kx_e = kx#np.sqrt(1+1j*1j*k*k+2*E/nu)#kp+1j*k##
	kx_h = -kx#np.sqrt(1+1j*1j*k*k-2*E/nu)#km-1j*k
		
	xiL_e = np.sqrt(2/nu)*(0+k*nu) #or -k*nu*sgn?
	xiL_h = np.sqrt(2/nu)*(0-k*nu) #or +k*nu*sgn?
	xiR_e = np.sqrt(2/nu)*(L+k*nu) #or -k*nu*sgn?
	xiR_h = np.sqrt(2/nu)*(L-k*nu) #or +k*nu*sgn?	
	
	pre_e1 = np.exp(-1j*nu*k+0.5*k**2*nu)
	pre_e2 = np.exp(1j*nu*k+0.5*k**2*nu)
	pre_h1 = np.exp(1j*nu*k+0.5*k**2*nu)
	pre_h2 = np.exp(-1j*nu*k+0.5*k**2*nu)
	
	D1_eL = pre_e1*Y1(a_e,xiL_e)
	D2_eL = pre_e2*Y2(a_e,xiL_e)
	D1_hL = pre_h1*Y1(a_h,xiL_h)
	D2_hL = pre_h2*Y2(a_h,xiL_h)
	
	D1_eR = pre_e1*Y1(a_e,xiR_e)
	D2_eR = pre_e2*Y2(a_e,xiR_e)
	D1_hR = pre_h1*Y1(a_h,xiR_h)
	D2_hR = pre_h2*Y2(a_h,xiR_h)
	
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
	
	D = det(matrix)
	return np.absolute(D)

def freeEnergy(phi,k,nu,delta,L,Z,kBT):
	n=10
	buf = 0.001*delta
	e0 = np.linspace(-delta+buf,delta-buf,n)
	E_array = []
	maxDiff = delta*10**(-6)
	for j in range(n):
		rootResult = opt.root(fun,e0[j],args=(k,nu,delta,L,Z,phi))
		diff = delta
		if len(E_array)>0:
			diff = rootResult - E_array[-1]
		if rootResult.success and diff>maxDiff:
			E_array.append(rootResult.x[0])
	result = 0
	for i in range(len(E_array)):
		result += np.log(2*np.cosh(0.5*E_array[i]/(delta*kBT)))
	return result

def dFreeEnergy(k,phi,nu,delta,L,Z,kBT):
	return misc.derivative(freeEnergy,phi,args=(k,nu,delta,L,Z,kBT),dx=0.001)
	
def totalCurrent(phi,nu,delta,L,Z,kBT):
	kmin = -1.
	kmin = 1.
	return integrate.quad(dFreeEnergy,kmin,kmax,args=(phi,nu,delta,L,Z,kBT)).y

def plotCurrent(nu,delta,L,Z,kBT,N):
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	I_array = np.zeros(phi_array.shape)
	for i in range(N):
		I_array[i] = totalCurrent(phi_arra[i],ni,delta,L,Z,kBT)
	plt.figure()
	plt.plot(phi_array,I_array,'.')
	plt.title('Total Current for nu ='+str(nu))
	plt.show()

nu = 200000.
delta = 200.
L = 106.7
kBT = delta
N = 100
plotCurrent(nu,delta,L,Z,kBT,N)
