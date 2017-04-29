import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy.optimize import brentq
from scipy.optimize import fsolve
from scipy.special import pbdv
from scipy.special import gamma
from scipy.linalg import det
from scipy.misc import factorial
import itertools
import scipy
import mpmath
from scipy import misc

def myexp(x):
	N = 10
	res = 0
	for i in range(N):
		res+=x**i/misc.factorial(i)
	return res

def myM2(a,xi):
	N = 40
	M = 0
	for i in range(N):
		M+=(np.sqrt(-1j*1j*a)*xi)**(2*i)/misc.factorial(2*i+1)
	return M
	
def myM1(a,xi):
	N = 40
	M = 0
	for i in range(N):
		M+=(np.sqrt(-1j*1j*a)*xi)**(2*i)/misc.factorial(2*i)
	return M
	

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
	#return np.exp(-0.25*xi**2)*M(0.5*a+0.25,0.5,0.5*xi**2)
	return myM1(a,xi)

def y2(a,xi):
	#return np.exp(-0.25*xi**2)*xi*M(0.5*a+0.75,1.5,0.5*xi**2)
	return myM2(a,xi)
	
def Y1(a,xi):
	return y1(a,xi)+np.sqrt(-1j*1j*a)*y2(a,xi)
	
def Y2(a,xi):
	return y1(a,xi)-np.sqrt(-1j*1j*a)*y2(a,xi)
	
def myFunc(ky,nu,a,xi):
	#return np.exp(-ky*np.sqrt(-1j*1j*2*a*nu)+0.5*ky**2*nu)*Y1(a,xi)
	return Y1(a,xi)
	
def limFunc(ky,nu,a,L,x,xi):
	k = np.sqrt(-2*a/nu)
	kx = k+1j*ky
	#return np.exp(1j*kx*(x+L/2))
	#return np.exp(-0.25*xi**2+np.sqrt(-1j*1j*a)*xi)
	return np.exp(np.sqrt(-1j*1j*a)*xi)
	

N = 100
L = 106.7
ky = 0.2
nu = np.logspace(0,10,N)
E = 2.
x = 0.1*L#np.linspace(-L/2,L/2,N)
myRes = np.zeros(nu.shape)
limRes = np.zeros(nu.shape)
dif = np.zeros(nu.shape)
xi = np.sqrt(2/nu)*((x+L/2)+ky*nu)
a = -(E+nu/2)

for i in range(N):
	myRes[i] = myFunc(ky,nu[i],a[i],xi[i]).imag
	limRes[i] = limFunc(ky,nu[i],a[i],L,x,xi[i]).imag
	dif[i] = np.absolute(myRes[i]-limRes[i])
	

plt.figure()
plt.plot(nu,myRes,label='myRes')
plt.plot(nu,limRes,label='limRes')
#plt.plot(nu,dif,label='dif')
#plt.semilogx(nu,dif,label='dif')
plt.xscale("log")
plt.yscale("log")
#plt.axis([100,1000,0,1])
plt.legend()
plt.show()