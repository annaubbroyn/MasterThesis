import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath
import math

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
	
	pre_e1 = 1#np.exp(-1j*nu*k+0.5*k**2*nu)
	pre_e2 = 1#np.exp(1j*nu*k+0.5*k**2*nu)
	pre_h1 = 1#np.exp(1j*nu*k+0.5*k**2*nu)
	pre_h2 = 1#np.exp(-1j*nu*k+0.5*k**2*nu)
	
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

def makePlotPhi(nu, delta, Z, k, L, N, n):
	print('makePlotPhi')
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
			
			#E_array[i] = opt.fsolve(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))[0]
		
		plt.plot(phi_array,E_array,'.b')
	plt.axis([phi_start,phi_end,-delta,delta])
	title = 'nu = '+str(nu)
	plt.title(title)
	plt.show()
	
def makePlotk(nu, delta, Z, phi, L, N, n):
	print('makePlotk')
	buf = 0.25
	k_start = -0.99
	k_end = 0.99
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

def testFunction(nu, delta, Z, k, L, N, n):
	print('testFunction')
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi = np.linspace(phi_start,phi_end,n)
	E_array = np.linspace(-delta,delta,N)
	Func_array = np.zeros(E_array.shape)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):		
			Func_array[i] = fun(E_array[i],k,nu,delta,L,Z,phi[j])
		plotLabel='phi = '+str(phi[j])
		plt.plot(E_array,Func_array,label=plotLabel)
	plt.legend()
	plt.show()
	
def plotAndSave(k_start, k_end, figCount, Z, k, L, N, n):
	print('plotAndSave')
	k_array = np.linspace(k_start,k_end,figCount)
	for k in k_array:
		delta = 2.
		nu = 2000.
		print('k:',k)
		print('nu:',nu)
		print('delta:',delta)
		buf = 0.1
		phi_start = -3*np.pi
		phi_end = 3*np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		E_array = np.zeros(phi_array.shape)
		e0 = np.linspace(-delta+buf,delta-buf,n)
		fig = plt.figure()
		for j in range(n):
			print('count:',j+1,'/',n)
			for i in range(N):	
				rootResult = opt.root(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))
				if rootResult.success:
					E_array[i] = rootResult.x[0]
				else:
					E_array[i] = 2*delta #this is a quick fix			
				#E_array[i] = opt.fsolve(fun,e0[j],args=(k,hw,delta,L,Z,phi_array[i]))[0]
			plt.plot(phi_array,E_array,'.b')
		plt.axis([phi_start,phi_end,-delta,delta])
		title = 'k = '+str(k)
		plt.title(title)
		#path = 'figures/042717/SNS_5_pcfd_delta_2_nu_200_root/'
		path = 'figures/042717/SNS_5_Only_y1_varying_k/delta_2_nu_2000/'
		name = 'k_'+str(k)+'n_'+str(n)+'N_'+str(N)
		fig.savefig(path+name+'.png')
	


k_start = 0.6
k_end = 0.95
figCount = 8
	
hw = 0.1#155789473684
N = 200
Z = 0
phi = 1.
L = 106.7

########
k0 = 0.
k2 = 0.2

delta0 = 200.
delta2 = 200.

n0 = 5
n2 = 5
########

k = k2
delta = 2.
n = 5
nu = 2000.

print('y1y2')
#print('k',k)
#print('delta',delta)
print('n',n)
#print('nu:',nu)
#makePlotPhi(nu,delta,Z,k,L,N,n)
#makePlotk(nu, delta, Z, phi, L, N, n)
#testFunction(nu, delta, Z, k, L, N, n)
plotAndSave(k_start, k_end, figCount, Z, k, L, N, n)

   
