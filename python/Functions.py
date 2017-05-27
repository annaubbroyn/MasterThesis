import numpy as np
from scipy import integrate
from scipy import special as sp
from scipy import misc
from scipy.linalg import det
from scipy import optimize as opt
import mpmath 
from scipy.interpolate import interp2d 

def y1(xi,a):
	res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
	factor = mpmath.exp(-0.25*xi**2)
	return complex(factor*res)
	
def y2(xi,a):
	res = xi*mpmath.hyp1f1(0.5*a+0.75,1.5,0.5*xi**2)
	factor = mpmath.exp(-0.25*xi**2)
	return complex(factor*res)
	
############################
print('Functions 1')

EF = 500.
Bmin = 1.
Bmax = 5.
L = 106.7
kymax = 0.2

amin = -(EF + 1)/Bmin
amax = -(EF - 1)/Bmin
anum = 50

xmin = -2*np.sqrt(EF/Bmin)*(kymax + Bmax*L/(2*EF))
xmax =  2*np.sqrt(EF/Bmin)*(kymax + Bmax*L/(2*EF))
#xmin = -np.sqrt(Bmax/EF)
#xmax = +np.sqrt(Bmax/EF) * (1 + Bmax*L/(2*EF))
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

print('Functions 2')

def Y1(x, a):
	return complex(y1_int_re(x, a)[0], y1_int_im(x,a)[0])
	
def Y2(x, a):
	return complex(y2_int_re(x, a)[0], y2_int_im(x,a)[0])

print('Functions 3')
	
############################
	

	
def fun(E_,phi,y,B,ky,Ef,L,Z,method):
	print('fun 1')
	
	
	E = E_[0]
	
	if E>1:
		E = 1.
	elif E<-1:
		E = -1.

	if abs(phi) == np.pi/2:
		phi +=0.01
	
		
	nu = 2*Ef/B
	beta = E/Ef
	alpha = L/nu
	
	a_e = -nu/2*(1+beta)
	a_h = -nu/2*(1-beta)	
		
	xiL_e = np.sqrt(2*nu)*(0+ky)
	xiL_h = np.sqrt(2*nu)*(0-ky)
	xiR_e = np.sqrt(2*nu)*(alpha+ky)
	xiR_h = np.sqrt(2*nu)*(alpha-ky)
	
	D1_eL = Y1(xiL_e,a_e)
	D2_eL = Y2(xiL_e,a_e)
	D1_hL = Y1(xiL_h,a_h)
	D2_hL = Y2(xiL_h,a_h)
	
	D1_eR = Y1(xiR_e,a_e)
	D2_eR = Y2(xiR_e,a_e)
	D1_hR = Y1(xiR_h,a_h)
	D2_hR = Y2(xiR_h,a_h)
	
	dD1_eL = np.sqrt(2/nu)*misc.derivative(Y1,xiL_e,args=(a_e,),dx=0.001)
	dD2_eL = np.sqrt(2/nu)*misc.derivative(Y2,xiL_e,args=(a_e,),dx=0.001)
	dD1_hL = np.sqrt(2/nu)*misc.derivative(Y1,xiL_h,args=(a_h,),dx=0.001)
	dD2_hL = np.sqrt(2/nu)*misc.derivative(Y2,xiL_h,args=(a_h,),dx=0.001)
	
	dD1_eR = np.sqrt(2/nu)*misc.derivative(Y1,xiR_e,args=(a_e,),dx=0.001)
	dD2_eR = np.sqrt(2/nu)*misc.derivative(Y2,xiR_e,args=(a_e,),dx=0.001)
	dD1_hR = np.sqrt(2/nu)*misc.derivative(Y1,xiR_h,args=(a_h,),dx=0.001)
	dD2_hR = np.sqrt(2/nu)*misc.derivative(Y2,xiR_h,args=(a_h,),dx=0.001)
	
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
					   [dD1_eL, dD2_eL, 0, 0, (ZL+1j)*qh*gamma_hL, (ZL-1j)*qe*gamma_eL, 0, 0],
					   [0, 0, dD1_hL, dD2_hL, (ZL+1j)*qh, (ZL-1j)*qe, 0, 0],
					   [D1_eR, D2_eR, 0, 0, 0, 0, gamma_hR, gamma_eR],
					   [0, 0, D1_hR, D2_hR, 0, 0, 1, 1],
					   [dD1_eR, dD2_eR, 0, 0, 0, 0, (-ZR-1j)*qh*gamma_hR, (-ZR+1j)*qe*gamma_eR],
					   [0, 0, dD1_hR, dD2_hR, 0, 0, (-ZR-1j)*qh, (-ZR+1j)*qe]])
	
	D = det(matrix)
	
	print('fun 2')
	
	return [D.real,D.imag]

def freeEnergy(phi,ky,y,B,Ef,L,Z,kBT,method):
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
		rootResult = opt.root(fun,e0[j],args=(phi,y,B,ky,Ef,L,Z,method))
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
	
	
def tabulate(Ef, B, ky, L, steps):
	E  = np.linspace(-1., 1., steps)

	nu = 2*Ef/B
	beta = E/Ef
	alpha = L/nu
	
	a_e = -nu/2*(1+beta)
	a_h = -nu/2*(1-beta)	
		
	xiL_e = np.sqrt(2*nu)*(0+ky)
	xiL_h = np.sqrt(2*nu)*(0-ky)
	xiR_e = np.sqrt(2*nu)*(alpha+ky)
	xiR_h = np.sqrt(2*nu)*(alpha-ky)
	
	D1_eL = interp(E, y1(xiL_e,a_e))
	D2_eL = interp(E, y2(xiL_e,a_e))
	D1_hL = interp(E, y1(xiL_h,a_h))
	D2_hL = interp(E, y2(xiL_h,a_h))
	
	D1_eR = interp(E, y1(xiR_e,a_e))
	D2_eR = interp(E, y2(xiR_e,a_e))
	D1_hR = interp(E, y1(xiR_h,a_h))
	D2_hR = interp(E, y2(xiR_h,a_h))
	
	return [D1_eL, D2_eL, D1_hL, D2_hL, D1_eR, D2_eR, D1_hR, D2_hR]

	
def dFreeEnergy(ky,phi,y,B,Ef,L,Z,kBT,method):
	print(ky)
	return  misc.derivative(freeEnergy,phi,args=(ky,y,B,Ef,L,Z,kBT,method),dx=0.001)
	
def currentDensity(y,phi,B,k_max,Ef,L,Z,kBT,method):
	if(k_max == 0):
		return dFreeEnergy(k_max,phi,y,B,Ef,L,Z,kBT,method)
	kyMin = -k_max
	kyMax = k_max
	return integrate.quad(dFreeEnergy,kyMin,kyMax,args=(phi,y,B,Ef,L,Z,kBT,method),limit=10)[0]

def totalCurrent(phi,B,k_max,Ef,L,W,Z,kBT,method,intLim):
	yMin = -W/2
	yMax = W/2
	return integrate.quad(currentDensity,yMin,yMax,args=(phi,B,k_max,Ef,L,Z,kBT,method),limit=intLim)[0]
