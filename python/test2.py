import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
from scipy import optimize as opt
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
import mpmath


def myExp(x,nu):
	xi = np.sqrt(2./nu)*np.absolute(x)
	a = (x/xi)**2
	M1 = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)#1+a*xi**2/2+(a+0.5)*(a+5/2)*xi**4/24#
	M2 = mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)#1+a*xi**2/6+(a+3/2)*(a+7/2)*xi**5/(5*24)#	
	return M1+x*M2

N = 500
x = np.linspace(0,100,N)
myRes = 1j*np.zeros(x.shape)
actRes = 1j*np.zeros(x.shape)
error = np.zeros(x.shape)
nu = np.logspace(3,7,N)
kFL = 106.7

for i in range(N):
	myRes[i] = myExp(1j*kFL,nu[i])
	actRes[i] = np.exp(1j*kFL)
	error[i] = np.absolute(myRes[i]-actRes[i])


plt.figure()
#plt.plot(nu,myRes.real,label='myRes')
#plt.plot(nu,actRes.real,label='actRes')
plt.plot(nu,error,label='error')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
