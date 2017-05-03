from Functions import fun
from Functions import totalCurrent
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')
import numpy as np
from scipy import optimize as opt

def makePlotEvsPhi(B,ky,Ef,L,Z,N,n,method):
	print('makePlotEvsPhi')
	print('method',method)
	print('ky',ky)
	print('Ef',Ef)
	print('B',B)
	print('N',N)
	
	de0 = 1/float(n+1)
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	E_array = np.zeros(phi_array.shape)
	e0 = np.linspace(-1+de0,1-de0,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):			
			rootResult = opt.root(fun,e0[j],args=(phi_array[i],B,ky,Ef,L,Z,method))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2
		plt.plot(phi_array,E_array,'.b')
	#plt.rc('text',usetex=True)
	plt.axis([phi_start,phi_end,-1,1])
	title = 'method = '+method + ', $\tilde{B}= '+str(B)+'$, $k_y/k_F = '+str(ky)+'$, $E_F/\Delta = '+str(Ef)+'$'
	plt.title(title)
	plt.show()

def makePlotEvsK(B,phi,Ef,L,Z,N,method):
	print('makePlotEvsK')
	print('method',method)
	print('phi',phi)
	print('Ef',Ef)
	print('B',B)
	print('n',n)
	print('N',N)
	
	ky_start = 0.
	ky_end = 0.5
	ky_array = np.linspace(ky_start, ky_end, N)
	E_array = np.zeros(ky_array.shape)
	plt.figure()
	for i in range(N):			
		rootResult = opt.root(fun,0.5,args=(phi,B,ky_array[i],Ef,L,Z,method))
		if rootResult.success:
			E_array[i] = rootResult.x[0]
		else:
			E_array[i] = 2
	plt.plot(ky_array,E_array,'.b')
	#plt.rc('text',usetex=True)
	plt.axis([ky_start,ky_end,0,1])
	title = 'method = '+method + ', $\tilde{B}= '+str(B)+'$, $phi = '+str(ky)+'$, $E_F/\Delta = '+str(Ef)+'$'
	plt.title(title)
	plt.show()

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
	
Z = 0
L = 106.7
N = 100
n = 2
Ef = 500.
ky = 0.0
B = 1.
phi = 0.
method = 'y1y2'

makePlotEvsPhi(B,ky,Ef,L,Z,100,n,method)
#makePlotEvsK(B,phi,Ef,L,Z,400,method)
   


   
