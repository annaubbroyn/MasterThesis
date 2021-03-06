import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath

##################################

def U(a,xi):
#Actually my y1
	return np.exp(-0.25*xi**2)*myM(0.5*a+0.25,0.5,0.5*xi**2)

def V(a,xi):
#Actually my y2
	return np.exp(-0.25*xi**2)*xi*myM(0.5*a+0.75,1.5,0.5*xi**2)

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
	
######################################

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
	
	dD1_eL = misc.derivative(dUFunction,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.000001)
	dD2_eL = misc.derivative(dVFunction,-L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.000001)
	dD1_hL = misc.derivative(dUFunction,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.000001)
	dD2_hL = misc.derivative(dVFunction,-L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.000001)
	
	dD1_eR = misc.derivative(dUFunction,L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.000001)
	dD2_eR = misc.derivative(dVFunction,L/2,args=(a_e,np.sqrt(2/nu),L/2+nu*k),dx=0.000001)
	dD1_hR = misc.derivative(dUFunction,L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.000001)
	dD2_hR = misc.derivative(dVFunction,L/2,args=(a_h,np.sqrt(2/nu),L/2-nu*k),dx=0.000001)
	
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
	return np.absolute(D)#D.real#D.real+D.imag
		
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
			"""
			rootResult = opt.root(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2*delta #this is a quick fix
			"""
			E_array[i] = opt.fsolve(fun,e0[j],args=(k,nu,delta,L,Z,phi_array[i]))[0]
		
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


		

hw = 1.#0.155789473684
N = 500
n = 3
nu = 40./hw
delta = 2./hw
Z = 0
phi = 1.
L = 106.7
k = 0.2

hw_start = 0.11
hw_end = 0.4
figCount = 20

print('Solve_SNS_2.py')
print('hw:',hw)
makePlotk(nu,delta,Z,phi,L,N,n)
#makePlotPhi(nu,delta,Z,k,L,N,n)
#testFunction(nu, delta, Z, k, L, N, n)
#plotAndSave(hw_start, hw_end, figCount, Z, k, L, N, n)
"""
a = 1
xi = 1
c1 = 1
c2 = 0
result = insertInEquation(a,xi,c1,c2)
print('result',result)


def UTemp(xi,a):
	return mpmath.pcfd(a,xi)
	
def VTemp(xi,a):
	return V(a,xi)
	
def insertInEquation(a,xi,c1,c2):
	y = c1*UTemp(xi,a)+c2*VTemp(xi,a)
	ddU = misc.derivative(UTemp,xi,args=(a,),n=2)
	ddV = misc.derivative(VTemp,xi,args=(a,),n=2)
	ddy = c1*ddU+c2*ddV
	return ddy-(xi*xi/4+a)*y
"""

   


