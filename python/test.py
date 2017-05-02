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


def y1(a,xi,k):
	nu = -2*a
	pre = 1#np.exp(-1j*nu*k+0.5*k**2*nu)
	res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
	factor = pre*np.exp(-0.25*xi**2)
	#This is a way to write factor*res without getting infs or nans
	return (float(factor.real*res.real-factor.imag*res.imag)+1j*float(factor.real*res.imag+factor.imag*res.real))
	
def y2(a,xi):
	res = xi*mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	factor = np.exp(-0.25*xi**2)
	return (float(factor.real*res.real-factor.imag*res.imag)+1j*float(factor.real*res.imag+factor.imag*res.real))
	
def Y1(a,xi,factor):
	#Low field:
	#return y1(a,xi)+np.sqrt(-1j*1j*a)*y2(a,xi)
	
	#Only y1 (fungerer for nu>100):
	pre = 1
	return pre*y1(a,xi,k)
	
	#y1 from pcfd
	#nu = -a*2
	#pre = np.exp(-1j*nu*k+0.5*k**2*nu)
	#Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	#dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	#C1 = (1-np.sin(np.pi*a))/(2*Ua0)
	#C2 = np.pi/(sp.gamma(0.5+a)*2*Ua0)
	#U = mpmath.pcfd(-a-1./2.,xi)
	#Um = mpmath.pcfd(-a-1./2.,xi)
	#V = (sp.gamma(1./2.+a)/np.pi)*(np.sin(np.pi*a)*U+Um)

	#High field:
	#Ua0 = np.sqrt(np.pi)/(2**(0.5*a+0.25)*sp.gamma(3./4.+0.5*a))
	#dUa0 = -np.sqrt(np.pi)/(2**(0.5*a-0.25)*sp.gamma(0.25+0.5*a))
	#return Ua0*y1(a,xi)+dUa0*y2(a,xi)
	
	#Pcfd
	#pre = 1
	#res = mpmath.pcfd(-a-1./2.,xi)
	#return float((pre*res).real)+1j*float((pre*res).imag)
	#return float(factor*mpmath.pcfd(-a-1./2.,xi))
	
	#Pcfu
	#pre = 1
	#res = mpmath.pcfu(a,xi)
	#print('res',res)
	#return float((pre*res).real)+1j*float((pre*res).imag)
	
def Y2(a,xi,factor):
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
	
	#Pcdf uten gamma
	#if np.isinf(sp.gamma(1./2.+a)):
	#	return np.inf
	#return np.sin(np.pi*a)*Y1(a,xi,factor)+Y1(a,-xi,factor)
	

	
nu = 2052.

a = -nu/2
k = 0.8
L = 106.7
xiL = np.sqrt(2/nu)*(L+k*nu)
xi0 = np.sqrt(2/nu)*k*nu
xi = xiL
factor = 10**(-0.43*nu)
print(Y1(a,xi,factor))
print('gamma(0.75+0.5a)',mpmath.gamma(3./4.+0.5*a))
