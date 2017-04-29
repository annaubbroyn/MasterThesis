import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath
import cmath

def U(a,x):
	#return float(sp.pbdv(-a-1./2.,x)[0])
    return float(mpmath.pcfd(-a-1./2.,x))

def V(a,x):
    if np.isinf(sp.gamma(1./2.+a)):
        return np.inf
    return (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U(a,x)+U(a,-x))

def mypcfd(a,z):
	return myU(-a-0.5,z)

def myU(a,z):
	U0a = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	dU0a = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	return np.exp(-0.25*z**2)*(U0a*M(0.5*a+0.25,0.5,0.5*z**2)+dU0a*z*M(0.5*a+3./4.,3./2.,0.5*z**2))

def myU2(a,b,z):
	return np.pi/np.sin(np.pi*b)*(M(a,b,z)/(sp.gamma(1+a-b)*sp.gamma(b))-z**(1-b)*M(1+a-b,2-b,z)/(sp.gamma(a)*sp.gamma(2-b)))

def M(a,b,z):
	return mpmath.hyp1f1(a,b,z)
	
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

def myy1(a,xi):
	return np.exp(-0.25*xi**2)*M(0.5*a+0.25,0.5,0.5*xi**2)

def myy2(a,xi):
	M = mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	return np.exp(-0.25*xi**2)*xi*(float(M.real)+1j*float(M.imag))
	#return np.exp(-0.25*xi**2)*xi*M(0.5*a+0.75,1.5,0.5*xi**2)
	
def dmyy1(xi,a):
	return myy1(a,xi)

def dmyy2(xi,a):
	return myy2(a,xi)

def y1(a,xi):
	Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	C1 = (1-np.sin(np.pi*a))/(2*Ua0)
	C2 = np.pi/(sp.gamma(0.5+a)*2*Ua0)
	return C1*U(a,xi)+C2*V(a,xi)

def y2(a,xi):
	Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	C1 = (1+np.sin(np.pi*a))/(2*dUa0)
	C2 = -np.pi/(sp.gamma(0.5+a)*2*dUa0)
	return C1*U(a,xi)+C2*V(a,xi)

def print1():
	N = 500	
	nu_start = 250
	nu_end = 350
	nu = np.linspace(nu_start,nu_end,N)
	mySol = np.zeros(nu.shape)
	libSol = np.zeros(nu.shape)
	E = 0
	L = 106.7
	x = 0
	a = -(E+nu/2)
	k = 0.2
	xi = np.sqrt(2/nu)*((x+L/2)+k*nu)
	print('xi_min',xi[0])
	print('xi_max',xi[-1])
	print('a_min',a[0])
	print('a_max',a[-1])

	for i in range(N):
		prefactor = np.exp(-(x+L/2)**2/(2*nu)+nu*k**2/2)
		mySol[i] = prefactor[i]*mypcfd(a[i],-xi[i])
		libSol[i] = prefactor[i]*mpmath.pcfd(a[i],-xi[i])
		
	plt.figure()
	plt.plot(nu,mySol,'.',label='mySol')
	plt.plot(nu,libSol,'.',label='libSol')
	#plt.axis([nu_start,nu_end,-5,5])	
	plt.legend()
	plt.show()


def print2():
	N = 500	
	#a = -18
	x = np.linspace(-3,3,N)
	a = -x**2
	mySol = np.zeros(x.shape)
	libSol = np.zeros(x.shape)
	for i in range(N):
		prefactor = 1
		mySol[i] = myU(a[i],x[i])
		libSol[i] = U(a[i],x[i])
		#mySol[i] = misc.derivative(myU,x[i],args=(x[i],))
		#libSol[i] = misc.derivative(U,x[i],args=(x[i],))
		
	plt.figure()
	plt.plot(x,mySol,'.',label='mySol')
	plt.plot(x,libSol,'.',label='libSol')
	#plt.axis([nu_start,nu_end,-5,5])	
	plt.legend()
	plt.show()

def print3():
	N = 1000	
	
	###
	#x = np.linspace(-5,5,N)
	#a = -18
	#mySol = np.zeros(x.shape)
	#libSol = np.zeros(x.shape)
	###
	nu = np.logspace(2,4,N)
	#nu = np.linspace(50,100,N)
	a = -nu/2
	L = 106.7
	xi = np.sqrt(2/nu)*L
	mySol = np.zeros(nu.shape)
	libSol = np.zeros(nu.shape)
	###
	
	for i in range(N):
		mySol[i] = myy2(a[i],xi[i])
		libSol[i] = y2(a[i],xi[i])
		#mySol[i] = misc.derivative(dmyy1,x[i],args=(a,))
		#libSol[i] = misc.derivative(dmyy1,x[i],args=(a,),dx=0.01)
		
		
	plt.figure()
	plt.plot(nu,mySol,'.',label='mySol')
	plt.plot(nu,libSol,'.',label='libSol')
	plt.xscale('log')
	#plt.yscale('log')
	#plt.axis([nu_start,nu_end,-5,5])	
	plt.legend()
	plt.show()

#print1()
#print2()
print3()