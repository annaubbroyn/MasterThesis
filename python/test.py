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
from scipy import integrate
from scipy.interpolate import interp2d 

def y1(xi,a):
	res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
	factor = mpmath.exp(-0.25*xi**2)
	return complex(factor*res)

def y2(xi,a):
	res = xi*mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	factor = mpmath.exp(-0.25*xi**2)
	return complex(factor*res)
	
EF = 500.
Bmin = 1.
Bmax = 8.
L = 106.7

amin = -(EF + 1)/Bmin
amax = -(EF - 1)/Bmin
anum = 50

xmin = -np.sqrt(Bmax/EF)
xmax = +np.sqrt(Bmax/EF) * (1 + Bmax*L/(2*EF))
xnum = 50

avec = np.linspace(amin, amax, anum)
xvec = np.linspace(xmin, xmax, xnum)
mxvec, mavec = np.meshgrid(xvec,avec)

y1_val = 1j*mavec
y2_val = 1j*mavec

for i, x in np.ndenumerate(mavec):
	y1_val[i] = y1(mxvec[i], mavec[i])
	y2_val[i] = y2(mxvec[i], mavec[i])
	
y1_int_re = interp2d(xvec, avec, y1_val.real, kind='cubic')
y1_int_im = interp2d(xvec, avec, y1_val.imag, kind='cubic')
y2_int_re = interp2d(xvec, avec, y2_val.real, kind='cubic')
y2_int_im = interp2d(xvec, avec, y2_val.imag, kind='cubic')

def Y1(x, a):
	return complex(y1_int_re(x, a)[0], y1_int_im(x,a)[0])
	
def Y2(x, a):
	return complex(y2_int_re(x, a)[0], y2_int_im(x,a)[0])


	
##########################################################################
	
y1real = np.zeros(100)
Y1real = np.zeros(100)
y2real = np.zeros(100)
Y2real = np.zeros(100)
dy1real = np.zeros(100)
dY1real = np.zeros(100)
dy2real = np.zeros(100)
dY2real = np.zeros(100)


"""
a = np.linspace(amin,amax,100)
x = (xmin+xmax)/2

for i in range(100):
	y1real[i] = y1(x,a[i]).real
	Y1real[i] = Y1(x,a[i]).real
	y2real[i] = y2(x,a[i]).real
	Y2real[i] = Y2(x,a[i]).real
	
plt.figure()
plt.plot(a,y1real)
plt.show()

plt.figure()
plt.plot(a,Y1real)
plt.show()

plt.figure()
plt.plot(a,y2real)
plt.show()

plt.figure()
plt.plot(a,Y2real)
plt.show()
"""

a = (amin+amax)/2
x = np.linspace(xmin,xmax,100)
for i in range(100):
	y1real[i] = y1(x[i],a).real
	Y1real[i] = Y1(x[i],a).real
	y2real[i] = y2(x[i],a).real
	Y2real[i] = Y2(x[i],a).real
	dy1real[i] = misc.derivative(y1,x[i],args=(a,),dx=0.001).real
	dY1real[i] = misc.derivative(Y1,x[i],args=(a,),dx=0.001).real
	dy2real[i] = misc.derivative(y2,x[i],args=(a,),dx=0.001).real
	dY2real[i] = misc.derivative(Y2,x[i],args=(a,),dx=0.001).real


#print('y1real',y1real)
#print('Y1real',Y1real)
#print('y2real',y2real)
#print('Y2real',Y2real)
"""
plt.figure()
plt.plot(x,y1real)
plt.show()

plt.figure()
plt.plot(x,Y1real)
plt.show()

plt.figure()
plt.plot(x,y2real)
plt.show()

plt.figure()
plt.plot(x,Y2real)
plt.show()
"""
plt.figure()
plt.plot(x,dy1real)
plt.show()

plt.figure()
plt.plot(x,dY1real)
plt.show()

plt.figure()
plt.plot(x,dy2real)
plt.show()

plt.figure()
plt.plot(x,dY2real)
plt.show()



	
