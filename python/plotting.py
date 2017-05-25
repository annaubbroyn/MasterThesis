from Functions import fun
from Functions import freeEnergy
from Functions import totalCurrent
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import time
import os
# nicer looking default plots
plt.style.use('bmh')
import numpy as np
from scipy import optimize as opt

############################################################
# E vs Phi

def plotAndSaveEvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,n,parabolic_method,variable):
	print('plotAndSaveEvsPhi')
	print('method',parabolic_method)
	if variable is not 'ky':
		print('ky',ky)
	if variable is not 'B':
		print('B',B)
	if variable is not 'Ef':
		print('Ef',Ef)
	print('N',N)
	print('n',n)
	x_array = np.linspace(start,end,figCount)
	for x in x_array:
		if variable is 'B':
			B = x
			print('B: ',B)
		elif variable is 'ky':
			ky = x
			print('ky: ',ky)
		elif variable is 'Ef':
			Ef = x
			print('Ef: ',Ef)
		
		
		de0 = 0.01
		e0 = np.linspace(-1+de0,1-de0,n)
		
		phi_start = -3*np.pi
		phi_end = 3*np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		
		E_array = np.zeros((2,N))
		time_array = np.zeros(N)
		success = np.zeros((n,N),dtype=bool)
		plotE_array = np.zeros((2,N))
		temp_E_array = np.zeros((n,N))
		success = np.zeros((n,N),dtype=bool)
		maxDiff = 10**(-6)
		
		for i in range(N):
			print('count:',i+1,'/',N)
			num_success = 0
			
			start = time.time()
			for j in range(n):
				rootResult = opt.root(fun,e0[j],args=(phi_array[i],y,B,ky,Ef,L,Z,parabolic_method))
				#rootResult = opt.root(fun,e0[j],args=(phi_array[i],B,ky,Ef,L,Z,parabolic_method),method = 'lm',options={'maxiter': 1000})
				temp_E_array[j][i] = float(rootResult.x[0])
				success[j][i] = rootResult.success
				if rootResult.success:
					num_success += 1
			if num_success == n:
				E_array[0][i] = temp_E_array[0][i]
				E_array[1][i] = temp_E_array[-1][i]
			elif num_success == 0:
				E_array[0][i] = -1.
				E_array[1][i] = 1.
			else:
				index = 0
				for j in range(n):
					if success[j][i]:
						E_array[0][i] = temp_E_array[j][i]
						index = j
						break
				for j in range(index,n):
					if (success[j][i] and np.absolute(E_array[0][i]-temp_E_array[j][i])>10**(-6)):
						E_array[1][i] = temp_E_array[j][i]
						break
					if j == n-1:
						if temp_E_array[0][i]>0.9999:
							E_array[1][i] = 1.
						elif E_array[0][i]<-0.9999:
							E_array[1][i] = -1.
						else:
							E_array[1][i] = E_array[0][i]
			E_array[0][i] = min(E_array[0][i],1.)
			E_array[0][i] = max(E_array[0][i],-1.)
			E_array[1][i] = min(E_array[1][i],1.)
			E_array[1][i] = max(E_array[1][i],-1.)
			
			end = time.time()
			time_array[i] = end-start;
		
		fig = plt.figure()
		plt.plot(phi_array,E_array[0],'.b')
		plt.plot(phi_array,E_array[1],'.b')
		plt.axis([phi_start,phi_end,-1-0.1,1+0.1])
		title = 'B = '+str(B) + ', Ef='+str(Ef) + ', ky='+str(ky)
		plt.title(title)		
		path = 'figures/050617/EvsPhi/Variable_'+variable+'/'
		folder = ""
		if variable is 'B':
			folder = 'Ef_%.1f_ky_%.3f/'%(Ef,ky)
		elif variable is 'Ef':
			folder = 'B_%.1f_ky_%.3f/'%(B,ky)
		if variable is 'ky':
			folder = 'Ef_%.1f_B_%.1f/'%(Ef,B)
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		name = 'method_%s_%s_%.2f_N_%d_n_%d' % (method,variable,x,N,n)
		name = name.replace('.','-')
		fig.savefig(path+name+'.png')
		plt.close(fig)
		
		fig = plt.figure()
		plt.plot(phi_array,time_array)
		title = 'Time vs phi for ABS energy'
		fig.savefig(path+name+'_time.png')
		plt.close(fig)
		
		
########################################################
# E vs k

def makePlotEvsKy(B,phi,Ef,L,Z,N,n,method):
	print('makePlotEvsKy')
	print('method',method)
	print('phi',phi)
	print('Ef',Ef)
	print('B',B)
	print('n',n)
	print('N',N)
	
	de0 = 0.05#1/float(n+1)
	ky_start = -1.25
	ky_end = 1.25
	ky_array = np.linspace(ky_start, ky_end, N)
	E_array = np.zeros(ky_array.shape)
	e0 = np.linspace(0+de0,1-de0,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):		
			rootResult = opt.root(fun,e0[j],args=(phi,B,ky_array[i],Ef,L,Z,method))
			if rootResult.success:
				E_array[i] = rootResult.x[0]
			else:
				E_array[i] = 2
		plt.plot(ky_array,E_array,'.b')
	plt.axis([ky_start,ky_end,0,1])	
	title = 'method = '+method + ', $\tilde{B}= '+str(B)+'$, $\phi = '+str(phi)+'$, $E_F/\Delta = '+str(Ef)+'$'
	plt.title(title)
	plt.show()



########################################################
# F vs Phi
	
def plotAndSaveFvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,method,variable):
	print('plotAndSaveFvsPhi')
	print('method',method)
	if variable is not 'ky':
		print('ky',ky)
	if variable is not 'B':
		print('B',B)
	if variable is not 'Ef':
		print('Ef',Ef)
	print('N',N)
	x_array = np.linspace(start,end,figCount)
	for x in x_array:
		if variable is 'B':
			B = x
			print('B: ',B)
		elif variable is 'ky':
			ky = x
			print('ky: ',ky)
		elif variable is 'Ef':
			Ef = x
			print('Ef: ',Ef)
		phi_start = -3*np.pi
		phi_end = 3*np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		F_array = np.zeros(phi_array.shape)
		time_array = np.zeros(phi_array.shape)
		for i in range(N):
			print('count:',i+1,'/',N)
			start = time.time()
			F_array[i] = freeEnergy(phi_array[i],ky,y,B,Ef,L,Z,kBT,method)
			end = time.time()
			time_array[i] = end-start
		
		#Plotting free energy
		fig = plt.figure()
		plt.plot(phi_array,F_array,'.')
		plt.axis([phi_start,phi_end,-2,-1])
		title = 'B = '+str(B) + ', Ef='+str(Ef) + ', ky='+str(ky)
		plt.title(title)
		path = 'figures/050617/FvsPhi/Variable_'+variable+'/'
		folder = ""
		if variable is 'B':
			folder = 'Ef_%.1f_ky_%.3f/'%(Ef,ky)
		elif variable is 'Ef':
			folder = 'B_%.1f_ky_%.3f/'%(B,ky)
		elif variable is 'ky':
			folder = 'Ef_%.1f_B_%.1f/'%(Ef,B)
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		name = 'method_%s_%s_%.2f_N_%d_y_%.2f_kmax_%.2f' % (method,variable,x,N,y,ky)
		name = name.replace('.','-')
		fig.savefig(path+name+'.png')
		plt.close(fig)
		
		#Plotting time
		fig = plt.figure()
		plt.plot(phi_array,time_array)
		title = 'Time vs phi for free energy'
		plt.title(title)
		fig.savefig(path+name+'_time.png')
		plt.close(fig)

########################################################
# Current vs Phi

def plotAndSaveCurrentvsPhi(start,end,figCount,k_max,B,Ef,L,W,Z,kBT,N,method,variable,intLim):
	print('plotAndSaveCurrentvsPhi')
	print('method',method)
	if variable is not 'B':
		print('B',B)
	if variable is not 'Ef':
		print('Ef',Ef)
	print('N',N)
	x_array = np.linspace(start,end,figCount)
	for x in x_array:
		if variable is 'B':
			B = x
			print('B: ',B)
		elif variable is 'ky':
			ky = x
			print('ky: ',ky)
		elif variable is 'Ef':
			Ef = x
			print('Ef: ',Ef)
		phi_start = -np.pi
		phi_end = np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		I_array = np.zeros(phi_array.shape)
		for i in range(N):
			print('count:',i+1,'/',N)
			print('      phi = ',phi_array[i])
			start = time.time()
			I_array[i] = totalCurrent(phi_array[i],B,k_max,Ef,L,W,Z,kBT,method,intLim)
			end = time.time()
			print('time spent: ',end-start)
			print(' ')
		print('test 1')
		fig = plt.figure()
		print('test 2')
		plt.plot(phi_array,I_array,'.')
		print('test 3')
		#plt.axis([phi_start,phi_end,-2,-1])
		title = 'B = '+str(B) + ', Ef='+str(Ef)
		plt.title(title)
		path = 'figures/050617/IvsPhi/withGauge/Variable_'+variable+'/'
		folder = ""
		if variable is 'B':
			folder = 'Ef_%.1f/k_max_%.2f/'% (Ef,k_max)
		elif variable is 'Ef':
			folder = 'B_%.1f/k_max_%.2f/'% (B,k_max)
		folder = folder.replace('.','-')
		path += folder
		print('test 4')
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		print('test 5')
		name = 'method_%s_%s_%.2f_N_%d_n_%d_W_%.2f_kmax_%.2f_limit_%d' % (method,variable,x,N,n,W,k_max,intLim)
		name = name.replace('.','-')
		print('test 6')
		fig.savefig(path+name+'.png')
		print('test 7')
		plt.close(fig)
		print('test 8')

		
def plotAndSaveCurrentvsB(B_start,B_end,k_max,Ef,L,W,Z,kBT,N,method,intLim):
	print('plotAndSaveCurrentvsB')
	print('method',method)
	print('Ef',Ef)
	print('N',N)
	print('k_max',k_max)
	B_array = np.linspace(B_start,B_end,N)
	phi = np.pi/2
	I_array = np.zeros(B_array.shape)
	for i in range(N):
		print('count:',i+1,'/',N)
		print('B: ',B_array[i])
		start = time.time()
		I_array[i] = totalCurrent(phi,B_array[i],k_max,Ef,L,W,Z,kBT,method,intLim)
		end = time.time()
		print('time spent: ',end-start)
		print(' ')
	fig = plt.figure()
	plt.plot(B_array,I_array,'.')
	#plt.axis([phi_start,phi_end,-2,-1])
	title = 'Current vs magnetic field'
	plt.title(title)
	path = 'figures/050617/IvsB/withGauge/'
	folder = 'Ef_%.1f/'%Ef
	folder = folder.replace('.','-')
	path += folder
	directory = os.path.dirname(path)
	if not os.path.exists(directory):
		os.makedirs(directory)
	name = 'method_%s_N_%d_Bstart_%.2f_Bend_%.2f_kMax_%.2f_W_%.2f_limit_%d' % (method,N,B_start,B_end,k_max,W,intLim)
	name = name.replace('.','-')
	fig.savefig(path+name+'.png')
	plt.close(fig)

########################################################
# Current vs B

	
def testFunction(B, Ef, ky, L, Z, N, n, method):
	print('testFunction')
	print('method',method)
	print('ky',ky)
	print('Ef',Ef)
	print('B',B)
	print('n',n)
	print('N',N)
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi = np.linspace(phi_start,phi_end,n)
	E_array = np.linspace(-1,1,N)
	Func_array = np.zeros(E_array.shape)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):		
			Func_array[i] = fun(E_array[i],phi[j],B,ky,Ef,L,Z,method)[0]
		plotLabel='phi = '+str(phi[j])
		plt.plot(E_array,Func_array,label=plotLabel)
	plt.legend()
	plt.show()
	
	
   	

intLim = 500
Z = 0
L = 106.7
W = L
N = 50
n = 4
Ef = 500.
ky = 0.
B = 2.3
phi = 1.#np.pi/2
kBT = 1.

method = 'y1y2'
#method = 'pcfuToy1y2'
#method = 'pcfu'

#variable = 'Ef'
#variable = 'ky'

variable = 'B'
start = 8.
end = 10.
k_max = 0.
figCount = 4
y = 4*Ef**2/(B**2*L)

plotAndSaveFvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,method,variable)
#plotAndSaveEvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,n,method,variable)
#makePlotEvsPhi(B,ky,Ef,L,Z,100,n,method)
#makePlotEvsKy(B,phi,Ef,L,Z,100,n,method)
#plotCurrentvsPhi(B,Ef,L,Z,kBT,N,method)
#testFunction(B, Ef, ky, L, Z, N, n, method)
#plotFvsPhi(B,Ef,ky,L,Z,kBT,N,method)
#plotAndSaveCurrentvsPhi(start,end,figCount,k_max,B,Ef,L,W,Z,kBT,N,method,variable,intLim)
#plotAndSaveCurrentvsB(start,end,k_max,Ef,L,W,Z,kBT,N,method,intLim)

#####################
#To remember
#In plot from papers: Ef = 10, B = 0.5
#In produced plot: Ef=500, B = 1.0 corresponds to nu=1000 and delta=1
   



