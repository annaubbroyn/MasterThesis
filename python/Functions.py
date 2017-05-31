import numpy as np
from scipy import integrate
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
from scipy import optimize as opt
import mpmath 
from scipy.interpolate import interp2d 
from copy import deepcopy as copy
import pickle
import os

def y1(xi,a):
	res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
	factor = mpmath.exp(-0.25*xi**2)
	return complex(factor*res)
	
def y2(xi,a):
	res = xi*mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	factor = mpmath.exp(-0.25*xi**2)
	return complex(factor*res)
	
############################

class parameters:
	def __init__(self,y = 0,ky = 0,phi=0,B=1,Bmin=0,Bmax=0,ky_max_int=0,ky_max_interp=0,anum=10,xnum=10,Ef=500,L=106.7,W=1067.,Z=0,kBT=1,interp = True):
		self.y = y
		self.ky = ky
		self.phi = phi
		self.B = B
		self.Ef = Ef
		self.L = L
		self.W = W
		self.Z = Z
		self.kBT = kBT
		self.anum = anum
		self.xnum = xnum
		self.Bmin = Bmin
		self.Bmax = Bmax
		self.ky_max_int = ky_max_int
		self.ky_max_interp = ky_max_interp
		self.interp = interp

class myObject:
	def __init__(self,param):
		if not param.interp:
			self.interp = False
			return
			
		self.interp = True
		if param.Bmin == 0:
			Bmin = param.B
		else:
			Bmin = param.Bmin
		if param.Bmax == 0:
			Bmax = param.B
		else:
			Bmax = param.Bmax
		if np.absolute(param.ky_max_interp) < np.absolute(param.ky):
			ky_max = np.absolute(param.ky)
		else:
			ky_max = param.ky_max_interp

		EF = param.Ef
		L = param.L
		anum = param.anum
		xnum = param.xnum

		amin = -(EF + 1)/Bmin
		amax = -(EF - 1)/Bmin

		xmin = -2*np.sqrt(EF/Bmin)*ky_max
		xmax =  2*np.sqrt(EF/Bmin)*(ky_max + Bmax*L/(2*EF))

		xmin -=0.1
		xmax +=0.1
		amin -=0.1
		amax +=0.1
		
		
		avec = np.linspace(amin, amax, anum)
		xvec = np.linspace(xmin, xmax, xnum)
		mxvec, mavec = np.meshgrid(xvec,avec)

		y1_val = 1j*mavec
		y2_val = 1j*mavec
		
		for i, x in np.ndenumerate(mavec):
			y1_val[i] = y1(mxvec[i], mavec[i])
			y2_val[i] = y2(mxvec[i], mavec[i])
		
		print('creating y1_int_re')
		self.y1_int_re = interp2d(xvec, avec, y1_val.real, kind='cubic')
		print('creating y1_int_im')
		self.y1_int_im = interp2d(xvec, avec, y1_val.imag, kind='cubic')
		print('creating y2_int_re')
		self.y2_int_re = interp2d(xvec, avec, y2_val.real, kind='cubic')
		print('creating y2_int_im')
		self.y2_int_im = interp2d(xvec, avec, y2_val.imag, kind='cubic')

class getMyObj:
	def __init__(self,param):
		filename = 'interpolation.bin'
		directory = os.path.dirname(filename)
		if not param.interp:
			self.interp = False
			return
		self.interp = True
		if os.path.isfile(filename):
			with open(filename,'rb') as f:
				self.y1_int_re, self.y1_int_im, self.y2_int_re, self.y2_int_im = pickle.load(f) 
		else:
			obj = myObject(param)
			self.y1_int_re = obj.y1_int_re
			self.y1_int_im = obj.y1_int_im
			self.y2_int_re = obj.y2_int_re
			self.y2_int_im = obj.y2_int_im

def Y1(x, a, obj):
	if obj.interp:
		return complex(obj.y1_int_re(x, a)[0], obj.y1_int_im(x,a)[0])
	return y1(x,a)
	
def Y2(x, a, obj):
	if obj.interp:
		return complex(obj.y2_int_re(x, a)[0], obj.y2_int_im(x,a)[0])
	return y2(x,a)
	
############################
	

	
def fun(E_,obj,param):
	
	E = E_[0]
	
	B = param.B
	Ef = param.Ef
	L = param.L
	Z = param.Z
	phi = param.phi
	ky = param.ky
	y = param.y
	
	if E>1:
		E = 1.
	elif E<-1:
		E = -1.

	"""
	if abs(phi) == np.pi/2:
		phi +=0.01
	"""
		
	nu = 2*Ef/B
	beta = E/Ef
	alpha = L/nu
	
	a_e = -nu/2*(1+beta)
	a_h = -nu/2*(1-beta)	
		
	xiL_e = np.sqrt(2*nu)*(0+ky)
	xiL_h = np.sqrt(2*nu)*(0-ky)
	xiR_e = np.sqrt(2*nu)*(alpha+ky)
	xiR_h = np.sqrt(2*nu)*(alpha-ky)
	
	D1_eL = Y1(xiL_e,a_e, obj)
	D2_eL = Y2(xiL_e,a_e, obj)
	D1_hL = Y1(xiL_h,a_h, obj)
	D2_hL = Y2(xiL_h,a_h, obj)
	
	D1_eR = Y1(xiR_e,a_e, obj)
	D2_eR = Y2(xiR_e,a_e, obj)
	D1_hR = Y1(xiR_h,a_h, obj)
	D2_hR = Y2(xiR_h,a_h, obj)
	
	dD1_eL = np.sqrt(2/nu)*misc.derivative(Y1,xiL_e,args=(a_e, obj),dx=0.001)
	dD2_eL = np.sqrt(2/nu)*misc.derivative(Y2,xiL_e,args=(a_e, obj),dx=0.001)
	dD1_hL = np.sqrt(2/nu)*misc.derivative(Y1,xiL_h,args=(a_h, obj),dx=0.001)
	dD2_hL = np.sqrt(2/nu)*misc.derivative(Y2,xiL_h,args=(a_h, obj),dx=0.001)
	
	dD1_eR = np.sqrt(2/nu)*misc.derivative(Y1,xiR_e,args=(a_e, obj),dx=0.001)
	dD2_eR = np.sqrt(2/nu)*misc.derivative(Y2,xiR_e,args=(a_e, obj),dx=0.001)
	dD1_hR = np.sqrt(2/nu)*misc.derivative(Y1,xiR_h,args=(a_h, obj),dx=0.001)
	dD2_hR = np.sqrt(2/nu)*misc.derivative(Y2,xiR_h,args=(a_h, obj),dx=0.001)
	
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
	
	D1_eL = D1_eL*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiL_e-ky))
	D2_eL = D2_eL*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiL_e-ky))
	D1_hL = D1_hL*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiL_h-ky))
	D2_hL = D2_hL*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiL_h-ky))
	D1_eR = D1_eR*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiR_e-ky))
	D2_eR = D2_eR*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiR_e-ky))
	D1_hR = D1_hR*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiR_h-ky))
	D2_hR = D2_hR*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiR_h-ky))
	dD1_eL = dD1_eL*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiL_e-ky))
	dD2_eL = dD2_eL*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiL_e-ky))
	dD1_hL = dD1_hL*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiL_h-ky))
	dD2_hL = dD2_hL*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiL_h-ky))
	dD1_eR = dD1_eR*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiR_e-ky))
	dD2_eR = dD2_eR*np.exp(1j*y*(0.5*np.sqrt(B/Ef)*xiR_e-ky))
	dD1_hR = dD1_hR*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiR_h-ky))
	dD2_hR = dD2_hR*np.exp(1j*y*(-0.5*np.sqrt(B/Ef)*xiR_h-ky))
	
	  
	
	matrix = np.array([[D1_eL, D2_eL, 0, 0, gamma_hL, gamma_eL, 0, 0],
					   [0, 0, D1_hL, D2_hL, 1, 1, 0, 0],
					   [dD1_eL, dD2_eL, 0, 0, (ZL*qh+1j/2*B/Ef*y+1j*qh)*gamma_hL, (ZL*qe+1j/2*B/Ef*y-1j*qe)*gamma_eL, 0, 0],
					   [0, 0, dD1_hL, dD2_hL, (ZL*qh-1j/2*B/Ef*y+1j*qh), (ZL*qe-1j/2*B/Ef*y-1j*qe), 0, 0],
					   [D1_eR, D2_eR, 0, 0, 0, 0, gamma_hR, gamma_eR],
					   [0, 0, D1_hR, D2_hR, 0, 0, 1, 1],
					   [dD1_eR, dD2_eR, 0, 0, 0, 0, (-ZR*qh+1j/2*B/Ef*y-1j*qh)*gamma_hR, (-ZR*qe+1j/2*B/Ef*y+1j*qe)*gamma_eR],
					   [0, 0, dD1_hR, dD2_hR, 0, 0, (-ZR*qh+1j/2*B/Ef*y-1j*qh), (-ZR*qe+1j/2*B/Ef*y+1j*qe)]])
	
	D = det(matrix)
	
	return [D.real,D.imag]

def freeEnergy(phi,obj,param):
	param.phi = phi
	n=4
	de0 = 0.01
	e0 = np.linspace(-1+de0,1-de0,n)
	temp_E_array = np.zeros(n)
	E_array = np.zeros(2)
	success = np.zeros(n,dtype=bool)
	maxDiff = 10**(-6)
	num_success = 0
	for j in range(n):
		rootResult = opt.root(fun,e0[j],args=(obj,param),method='lm',options={'xtol': 1e-08,})
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
		index = 0
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
				if temp_E_array[0]>1.:
					E_array[1] = 1.
				elif temp_E_array[0]<-1.:
					E_array[1] = -1.
				else:
					E_array[1] = E_array[0]
	E_array[0] = min(E_array[0],1.)
	E_array[0] = max(E_array[0],-1.)
	E_array[1] = min(E_array[1],1.)
	E_array[1] = max(E_array[1],-1.)
			
	result = 0
	for i in range(2):
		result += -np.log(2*np.cosh(0.5*E_array[i]/param.kBT))
	return result					

	
def dFreeEnergy(ky,param):
	param_copy = copy(param)
	param_copy.ky = ky
	phi = param_copy.phi
	obj=getMyObj(param_copy)
	return  misc.derivative(freeEnergy,phi,args=(obj,param_copy),dx=0.1)
	
def currentDensity(y,param):
	param_copy = copy(param)
	param_copy.y = y
	k_max = param_copy.ky_max_int
	#print('phi in totalCurrent',param_copy.phi)
	if(k_max == 0):
		return dFreeEnergy(0.,param_copy)
	kyMin = -k_max
	kyMax = k_max
	return integrate.quad(dFreeEnergy,kyMin,kyMax,args=(param_copy),limit=50)[0]

def totalCurrent(param):
	copy_param = copy(param)
	if param.W == 0:
		return currentDensity(copy_param.y,copy_param)
	yMin = -copy_param.W/2
	yMax = copy_param.W/2
	print('phi in totalCurrent',copy_param.phi)
	return integrate.quad(currentDensity,yMin,yMax,args=(copy_param,),limit=50)[0]


	
	
	
	
	