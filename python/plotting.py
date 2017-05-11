from Functions import fun
from Functions import freeEnergy
from Functions import totalCurrent
from matplotlib import pyplot as plt
import time
import os
# nicer looking default plots
plt.style.use('bmh')
import numpy as np
from scipy import optimize as opt

############################################################
# E vs Phi

def makePlotEvsPhi(B,ky,Ef,L,Z,N,n,method):
	print('makePlotEvsPhi')
	print('method',method)
	print('ky',ky)
	print('Ef',Ef)
	print('B',B)
	print('N',N)
	print('n',n)
	
	de0 = 0.05#2/float(n+1)
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	E_array = np.zeros(phi_array.shape)
	time_array = np.zeros(phi_array.shape)
	e0 = np.linspace(-1+de0,1-de0,n)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):
			start = time.time()
			rootResult = opt.root(fun,e0[j],args=(phi_array[i],B,ky,Ef,L,Z,method))
			end = time.time()
			time_array[i] = end-start
			#if rootResult.success:
			E_array[i] = rootResult.x[0]
			#else:
			#	E_array[i] = 2
		plt.plot(phi_array,E_array,'.b')
	#plt.rc('text',usetex=True)
	plt.axis([phi_start,phi_end,-1,1])
	title = 'method = '+method + ', $\tilde{B}= '+str(B)+'$, $k_y/k_F = '+str(ky)+'$, $E_F/\Delta = '+str(Ef)+'$'
	plt.title(title)
	plt.show()
	plt.figure()
	plt.plot(phi_array,time_array)
	plt.show()
	
def plotAndSaveEvsPhi(start,end,figCount,B,ky,Ef,L,Z,N,n,method,variable):
	print('plotAndSaveEvsPhi')
	print('method',method)
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
		de0 = 0.#2/float(n+1)
		phi_start = -3*np.pi
		phi_end = 3*np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		E_array = np.zeros(phi_array.shape)
		temp_E_array = np.zeros(phi_array.shape)
		#plot_temp = np.zeros(phi_array.shape)
		dE_array = np.zeros(phi_array.shape)
		e0 = np.linspace(-1+de0,1-de0,n)
		#e0 = [-1+de0]
		fig = plt.figure()
		for j in range(n):
			print('count:',j+1,'/',n)
			for i in range(N):
				if i == 0:
					e0_now = e0[j]
					rootResult = opt.root(fun,e0_now,args=(phi_array[i],B,ky,Ef,L,Z,method))
					if e0_now > 0.:
						while not rootResult.success and e0_now>=-1.:
							e0_now -= 0.2
							rootResult = opt.root(fun,e0_now,args=(phi_array[i],B,ky,Ef,L,Z,method))
					if e0_now < 0.:
						while not rootResult.success and e0_now<=1.:
							e0_now += 0.2
							rootResult = opt.root(fun,e0_now,args=(phi_array[i],B,ky,Ef,L,Z,method))
				else:
					e0_now = temp_E_array[i-1]+dE_array[i-1]
					rootResult = opt.root(fun,e0_now,args=(phi_array[i],B,ky,Ef,L,Z,method))
				
				#if i is 0 or rootResult.success:
				temp_E_array[i] = rootResult.x[0]
				#else:
				#	temp_E_array[i] = temp_E_array[i-1]
				
				aboveZero = False
				dphi=(phi_end-phi_start)/N
				index = int((phi_array[i]-np.pi-phi_start)/dphi) + 5;
				if temp_E_array[index] > 0.:
					aboveZero = True

				if i == 0:
					#print('a')
					dE_array[i] = 0.
					#temp_E_array[i] = min(0.9998,temp_E_array[i])
					#temp_E_array[i] = max(-0.9998,temp_E_array[i])
				elif temp_E_array[i]>=0.99999:
					#print('b')
					#print('de_array[i-1] is ',dE_array[i-1])
					if dE_array[i-1]>0.0 and aboveZero:
						#print('is in if')
						temp_E_array[i] = -1.
						dE_array[i] = 0.0001
					else:
						#print('is in else')
						#while not rootResult.success and e0_now>-1. and temp_E_array[i-2] > 0.:
						#	e0_now -= 0.1
						#	rootResult = opt.root(fun,e0_now,args=(phi_array[i],B,ky,Ef,L,Z,method))
						#if rootResult.success:
						#	temp_E_array[i] = rootResult.x[0]
						#	dE_array[i]=rootResult.x[0]-temp_E_array[i-1]
						#else:
						temp_E_array[i] = 1.
						dE_array[i] = -0.0001
				elif temp_E_array[i]<=-0.99999:
					#print('c')
					if dE_array[i-1]<0.0 and not aboveZero:
						temp_E_array[i] = 1.
						dE_array[i] = -0.0001
					else:
						#while not rootResult.success and e0_now<1. and temp_E_array[i-2] < 0.:
						#	e0_now += 0.1
						#	rootResult = opt.root(fun,e0_now,args=(phi_array[i],B,ky,Ef,L,Z,method))
						#if rootResult.success:
						#	temp_E_array[i] = rootResult.x[0]
						#	dE_array[i]=rootResult.x[0]-temp_E_array[i-1]
						#else:
						temp_E_array[i] = -1.
						dE_array[i] = 0.0001
				elif rootResult.success:
					dE_array[i] = rootResult.x[0]-temp_E_array[i-1]
				else:
					#print('d')
					temp_E_array[i] = temp_E_array[i-1]
					dE_array[i] = 0.
				
				dE_array[i] = min(dE_array[i],0.05)
				dE_array[i] = max(dE_array[i],-0.05)
				
				if rootResult.success:
					E_array[i] = rootResult.x[0]
				else:
					E_array[i] = 2 #this is a quick fix	
					#print('at phi =',phi_array[i], ' rootResult is',rootResult.x[0], ', temp_E_array is ', temp_E_array[i], 'and dE_array[i] is', dE_array[i])
				"""
				if j == 0:
					plot_temp[i] = temp_E_array[i]
				
				if phi_array[i]>-5:
				
				#if j == 1:				
					print('E_array[i]',E_array[i])
					print('temp_E_array[i]',temp_E_array[i])
					print('dE_array[i]',dE_array[i])
					print('rootResult.x[0]',rootResult.x[0])
					print('success?',rootResult.success)
					print(' ')
				"""
			plt.plot(phi_array,E_array,'.b')
			#print(temp_E_array)
		plt.axis([phi_start,phi_end,-1,1])
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
		
		#plt.figure()
		#plt.plot(phi_array,plot_temp)
		#plt.show()
		
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

def plotFvsPhi(B,Ef,ky,L,Z,kBT,N,method):
	print('makePlotFvsPhi')
	print('method',method)
	print('Ef',Ef)
	print('B',B)
	print('ky',ky)
	print('kBT',kBT)
	print('N',N)
	
	phi_start = -3*np.pi
	phi_end = 3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	F_array = np.zeros(phi_array.shape)
	time_array = np.zeros(phi_array.shape)
	for i in range(N):
		start = time.time()
		F_array[i] = freeEnergy(phi_array[i],ky,B,Ef,L,Z,kBT,method)
		end = time.time()
		time_array[i] = end-start
	plt.figure()
	plt.plot(phi_array,F_array,'.')
	title = 'Total Current with method = '+method + ', $\tilde{B}= '+str(B)+'$, $k_y/k_F = '+str(ky)+'$, $E_F/\Delta = '+str(Ef)+'$'
	#plt.axis([phi_start,phi_end,-2,-1])
	plt.show()
	plt.figure()
	plt.plot(phi_array,time_array)
	title = 'Calculation time'
	plt.show()
	
def plotAndSaveFvsPhi(start,end,figCount,B,ky,Ef,L,Z,N,method,variable):
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
		for i in range(N):
			F_array[i] = freeEnergy(phi_array[i],ky,B,Ef,L,Z,kBT,method)
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
		name = 'method_%s_%s_%.2f_N_%d' % (method,variable,x,N)
		name = name.replace('.','-')
		fig.savefig(path+name+'.png')
		plt.close(fig)



########################################################
# Current vs Phi
def plotCurrentvsPhi(B,Ef,L,Z,kBT,N,method):
	print('makePlotCurrentvsPhi')
	print('method',method)
	print('Ef',Ef)
	print('B',B)
	print('kBT',kBT)
	print('N',N)
	phi_start = -3*np.pi
	phi_end = -7.1#3*np.pi
	phi_array = np.linspace(phi_start, phi_end, N)
	I_array = np.zeros(phi_array.shape)
	for i in range(N):
		print('count: ',i+1,'/',N)
		I_array[i] = totalCurrent(phi_array[i],B,Ef,L,Z,kBT,method)
	plt.figure()
	plt.plot(phi_array,I_array,'.')
	title = 'Total Current with method = '+method + ', $\tilde{B}= '+str(B)+'$, $k_y/k_F = '+str(ky)+'$, $E_F/\Delta = '+str(Ef)+'$'
	plt.show()
	
def plotAndSaveCurrentvsPhi(start,end,figCount,B,Ef,L,Z,kBT,N,method,variable):
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
		phi_start = -3*np.pi
		phi_end = 3*np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		I_array = np.zeros(phi_array.shape)
		for i in range(N):
			print('count:',i+1,'/',N)
			I_array[i] = totalCurrent(phi_array[i],B,Ef,L,Z,kBT,method)
		fig = plt.figure()
		plt.plot(phi_array,I_array,'.')
		plt.axis([phi_start,phi_end,-2,-1])
		title = 'B = '+str(B) + ', Ef='+str(Ef)
		plt.title(title)
		path = 'figures/050617/IvsPhi/Variable_'+variable+'/'
		folder = ""
		if variable is 'B':
			folder = 'Ef_%.1f/'%Ef
		elif variable is 'Ef':
			folder = 'B_%.1f/'%B
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		name = 'method_%s_%s_%.2f_N_%d_n_%d' % (method,variable,x,N,n)
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
	
	
   	

	
Z = 0
L = 106.7
N = 100
n = 2
Ef = 500.
ky = 0.2
B = 2.3
phi = 1.#np.pi/2
kBT = 1.

method = 'y1y2'
#method = 'pcfuToy1y2'
#method = 'pcfu'

#variable = 'Ef'
#variable = 'ky'

variable = 'B'
start = 0.1
end = 10.
figCount = 10

#plotAndSaveFvsPhi(start,end,figCount,B,ky,Ef,L,Z,N,method,variable)
plotAndSaveEvsPhi(start,end,figCount,B,ky,Ef,L,Z,N,n,method,variable)
#makePlotEvsPhi(B,ky,Ef,L,Z,100,n,method)
#makePlotEvsKy(B,phi,Ef,L,Z,100,n,method)
#plotCurrentvsPhi(B,Ef,L,Z,kBT,N,method)
#testFunction(B, Ef, ky, L, Z, N, n, method)
#plotFvsPhi(B,Ef,ky,L,Z,kBT,N,method)
#plotAndSaveCurrentvsPhi(start,end,figCount,B,Ef,L,Z,kBT,N,method,variable)

#####################
#To remember
#In plot from papers: Ef = 10, B = 0.5
#In produced plot: Ef=500, B = 1.0 corresponds to nu=1000 and delta=1
   


   
