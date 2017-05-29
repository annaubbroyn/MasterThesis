import Functions
import matplotlib
matplotlib.use('Agg')
from Functions import fun
from Functions import freeEnergy
from Functions import totalCurrent
from Functions import totalCurrent2
from Functions import currentDensity
from Functions import dFreeEnergy
from matplotlib import pyplot as plt
from Functions import parameters
from Functions import myObject
import time
import os
# nicer looking default plots
plt.style.use('bmh')
import numpy as np
from scipy import optimize as opt
from copy import deepcopy as copy

############################################################
# E vs Phi

def plotAndSaveEvsPhi(variable,start,end,figCount,xVariable,startXval,endXval,param,N,n):
	print('plotAndSaveEvs',xVariable)
	if variable is not 'ky' and xVariable is not 'ky':
		print('ky',param.ky)
	if variable is not 'B' and xVariable is not 'B':
		print('B',param.B)
	if variable is not 'Ef' and xVariable is not 'Ef':
		print('Ef',param.Ef)
	if variable is not 'y' and xVariable is not 'y':
		print('y',param.y)
	if variable is not 'phi' and xVariable is not 'phi':
		print('phi',param.phi)
	print('N',N)
	print('n',n)
	x_array = np.linspace(start,end,figCount)
	for x in x_array:
		if variable is 'B':
			param.B = x
			print('B: ',x)
		elif variable is 'ky':
			param.ky = x
			print('ky: ',x)
		elif variable is 'Ef':
			param.Ef = x
			print('Ef: ',x)
		elif variable is 'y':
			param.y = x
			print('y: ',x)
		elif variable is 'phi':
			param.phi = x
			print('phi: ',x)
		elif variable is 'Bmax':
			param.Bmax = x
			print('Bmax: ',x)
		elif variable is 'Bmin':
			param.Bmin = x
			print('Bmin: ',x)
		elif variable is 'kyInterp':
			param.kyInterp = x
			print('kyInterp: ',x)
		elif variable is 'anum':
			param.anum = x
			print('anum: ',x)
		elif variable is 'xnum':
			param.xnum = x
			print('xnum: ',x)
		
		path = 'figures/052717/Evs'+xVariable+'/figVariable_'+variable+'/'
		folder = ""
		title = 'E vs '+xVariable + ', ' + variable + ' = ' + str(x)
		
		if variable is not 'B' and xVariable is not 'B':
			folder+='B_%.1f_'%param.B
			title += ' B = %.1f' %param.B
		if variable is not 'Ef' and xVariable is not 'Ef':
			folder+='Ef_%.1f_'%param.Ef
			title += ' Ef = %.1f' %param.Ef
		if variable is not 'ky' and xVariable is not 'ky':
			folder+='ky_%.3f_'%param.ky
			title += ' ky = %.3f' %param.ky
		if variable is not 'y' and xVariable is not 'y':
			folder+='y_%.1f_'%param.y		
			title += ' y = %.1f' %param.y			
		if variable is not 'phi' and xVariable is not 'phi':
			folder+='phi_%.1f_'%param.phi			
			title += ' phi = %.1f' %param.phi
		if variable is not 'Bmax' and xVariable is not 'Bmax':
			folder+='Bmax_%.1f_'%param.Bmax			
			title += ' Bmax = %.1f' %param.Bmax
		if variable is not 'Bmin' and xVariable is not 'Bmin':
			folder+='Bmin_%.1f_'%param.Bmin			
			title += ' Bmin = %.1f' %param.Bmin
		if variable is not 'kyInterp' and xVariable is not 'kyInterp':
			folder+='kyInterp_%.1f_'%param.ky_max_interp		
			title += 'kyInterp = %.1f' %param.ky_max_interp		
		if variable is not 'anum' and xVariable is not 'anum':
			folder+='anum_%.1f_'%param.anum		
			title += 'anum = %.1f' %param.anum
		if variable is not 'xnum' and xVariable is not 'xnum':
			folder+='xnum_%.1f_'%param.xnum		
			title += 'xnum = %.1f' %param.xnum
		
		folder += '/'
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
			
		prefactor = '%.2f' % (100-x)
		prefactor = prefactor.replace('.','-')
		name = '%s_%.2f_N_%d_n_%d_from_%.2f_to_%.2f' % (variable,x,N,n,startXval,endXval)
		name = name.replace('.','-')
		if param.interp:
			name+='_interp'
		name = prefactor+name
		
		
		de0 = 0.01
		e0 = np.linspace(-1+de0,1-de0,n)
		
		xVal_array = np.linspace(startXval,endXval, N)
		
		E_array = np.zeros((2,N))
		time_array = np.zeros(N)
		success = np.zeros((n,N),dtype=bool)
		plotE_array = np.zeros((2,N))
		temp_E_array = np.zeros((n,N))
		success = np.zeros((n,N),dtype=bool)
		maxDiff = 10**(-6)
		
		if xVariable is not 'B':
			start_interp_time = time.time()
			obj = myObject(param)
			end_interp_time = time.time()
			time_interp = end_interp_time - start_interp_time
		
		for i in range(N):
			print('count:',i+1,'/',N)
			num_success = 0
			
			if xVariable is 'y':
				param.y = xVal_array[i]
			elif xVariable is 'phi':
				param.phi=xVal_array[i]
			elif xVariable is 'B':
				param.B = xVal_array[i]
				start_interp_time = time.time()
				obj = myObject(param)
				end_interp_time = time.time()
				time_interp = end_interp_time - start_interp_time
			else:
				print('xVariable is unvalid')
				return
			

			
			for j in range(n):
			
				start = time.time()			
				rootResult = opt.root(fun,e0[j],args=(obj,param))
				end = time.time()
				time_array[i] += end-start;
			
				temp_E_array[j][i] = float(rootResult.x[0])
				success[j][i] = rootResult.success
				if rootResult.success:
					num_success += 1

			time_array[i]=time_array[i]/n

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
			
		
		fig = plt.figure()
		plt.plot(xVal_array,E_array[0],'.b')
		plt.plot(xVal_array,E_array[1],'.b')
		plt.axis([startXval,endXval,-1-0.1,1+0.1])
		plt.title(title)
		fig.savefig(path+name+'.png')
		plt.close(fig)
		
		
		fig = plt.figure()
		plt.plot(xVal_array,time_array)
		title = 'Time vs %s for ABS energy. Interpolation time %.2f' %(xVariable,time_interp)
		plt.title(title)
		fig.savefig(path+'time_'+name+'.png')
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
def plotAndSaveFvsPhi(variable,start,end,figCount,xVariable,startXval,endXval,param,N):
#def plotAndSaveFvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,method,variable):
	print('plotAndSaveFvs',xVariable)
	if variable is not 'ky' and xVariable is not 'ky':
		print('ky',param.ky)
	if variable is not 'B' and xVariable is not 'B':
		print('B',param.B)
	if variable is not 'Ef' and xVariable is not 'Ef':
		print('Ef',param.Ef)
	if variable is not 'y' and xVariable is not 'y':
		print('y',param.y)
	if variable is not 'phi' and xVariable is not 'phi':
		print('phi',param.phi)
	print('N',N)
	x_array = np.linspace(start,end,figCount)
	
	
	for x in x_array:
		if variable is 'B':
			param.B = x
			print('B: ',x)
		elif variable is 'ky':
			param.ky = x
			print('ky: ',x)
		elif variable is 'Ef':
			param.Ef = x
			print('Ef: ',x)
		elif variable is 'y':
			param.y = x
			print('y: ',x)
		elif variable is 'phi':
			param.phi = x
			print('phi: ',x)
		elif variable is 'Bmax':
			param.Bmax = x
			print('Bmax: ',x)
		elif variable is 'Bmin':
			param.Bmin = x
			print('Bmin: ',x)
		elif variable is 'kyInterp':
			param.kyInterp = x
			print('kyInterp: ',x)
		elif variable is 'anum':
			param.anum = x
			print('anum: ',x)
		elif variable is 'xnum':
			param.xnum = x
			print('xnum: ',x)
			
		path = 'figures/052717/Fvs'+xVariable+'/figVariable_'+variable+'/'
		folder = ""
		title = 'F vs '+xVariable + ', ' + variable + ' = ' + str(x)
		
		if variable is not 'B' and xVariable is not 'B':
			folder+='B_%.1f_'%param.B
			title += ' B = %.1f' %param.B
		if variable is not 'Ef' and xVariable is not 'Ef':
			folder+='Ef_%.1f_'%param.Ef
			title += ' Ef = %.1f' %param.Ef
		if variable is not 'ky' and xVariable is not 'ky':
			folder+='ky_%.3f_'%param.ky
			title += ' ky = %.3f' %param.ky
		if variable is not 'y' and xVariable is not 'y':
			folder+='y_%.1f_'%param.y		
			title += ' y = %.1f' %param.y			
		if variable is not 'phi' and xVariable is not 'phi':
			folder+='phi_%.1f_'%param.phi			
			title += ' phi = %.1f' %param.phi
		if variable is not 'Bmax' and xVariable is not 'Bmax':
			folder+='Bmax_%.1f_'%param.Bmax			
			title += ' Bmax = %.1f' %param.Bmax
		if variable is not 'Bmin' and xVariable is not 'Bmin':
			folder+='Bmin_%.1f_'%param.Bmin			
			title += ' Bmin = %.1f' %param.Bmin
		if variable is not 'kyInterp' and xVariable is not 'kyInterp':
			folder+='kyInterp_%.1f_'%param.ky_max_interp		
			title += 'kyInterp = %.1f' %param.ky_max_interp		
		if variable is not 'anum' and xVariable is not 'anum':
			folder+='anum_%.1f_'%param.anum		
			title += 'anum = %.1f' %param.anum
		if variable is not 'xnum' and xVariable is not 'xnum':
			folder+='xnum_%.1f_'%param.xnum		
			title += 'xnum = %.1f' %param.xnum
		
		
		folder += '/'
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		prefactor = '%.2f' % (100-x)
		prefactor = prefactor.replace('.','-')
		name = '%s_%.2f_N_%d_n_%d_from_%.2f_to_%.2f' % (variable,x,N,n,startXval,endXval)
		name = name.replace('.','-')
		if param.interp:
			name+='_interp'
		name = prefactor + name
		
		
		if xVariable is not 'B':
			start_interp_time = time.time()
			obj = myObject(param)
			end_interp_time = time.time()
			time_interp = end_interp_time - start_interp_time
		
		xVal_array = np.linspace(startXval,endXval, N)
		F_array = np.zeros(xVal_array.shape)
		time_array = np.zeros(xVal_array.shape)
		for i in range(N):
			if xVariable is 'y':
				param.y = xVal_array[i]
			elif xVariable is 'phi':
				param.phi=xVal_array[i]
			elif xVariable is 'ky':
				param.ky=xVal_array[i]
			elif xVariable is 'B':
				param.B = xVal_array[i]
				start_interp_time = time.time()
				obj = myObject(param)
				end_interp_time = time.time()
				time_interp = end_interp_time - start_interp_time
			else:
				print('xVariable is unvalid')
				return -1
		
			print('count:',i+1,'/',N)
			start = time.time()
			F_array[i] = freeEnergy(param.phi,obj,param)
			end = time.time()
			time_array[i] = end-start
		
		fig = plt.figure()
		plt.plot(xVal_array,F_array,'.')
		#plt.axis([startXval,endXval,-2,-1])
		plt.title(title)
		fig.savefig(path+name+'.png')
		plt.close(fig)
		
		fig = plt.figure()
		plt.plot(xVal_array,time_array)
		title = 'Time vs %s for Free energy. Interpolation time %.2f' %(xVariable,time_interp)
		plt.title(title)
		fig.savefig(path+'time_'+name+'.png')
		plt.close(fig)

		

##########################################################
def plotAndSaveCurrentDensityvsy(variable,start,end,figCount,starty,endy,param,N):
#def plotAndSaveFvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,method,variable):
	print('plotAndSaveCurrentDensityvsy')
	if variable is not 'ky':
		print('ky',param.ky)
	if variable is not 'B':
		print('B',param.B)
	if variable is not 'Ef':
		print('Ef',param.Ef)
	if variable is not 'phi':
		print('phi',param.phi)
	print('N',N)	
	
	x_array = np.linspace(start,end,figCount)
	
	
	for x in x_array:
		if variable is 'B':
			param.B = x
			print('B: ',x)
		elif variable is 'ky':
			param.ky = x
			print('ky: ',x)
		elif variable is 'Ef':
			param.Ef = x
			print('Ef: ',x)
		elif variable is 'y':
			param.y = x
			print('y: ',x)
		elif variable is 'phi':
			param.phi = x
			print('phi: ',x)
		elif variable is 'Bmax':
			param.Bmax = x
			print('Bmax: ',x)
		elif variable is 'Bmin':
			param.Bmin = x
			print('Bmin: ',x)
		elif variable is 'kyInterp':
			param.kyInterp = x
			print('kyInterp: ',x)
		elif variable is 'anum':
			param.anum = x
			print('anum: ',x)
		elif variable is 'xnum':
			param.xnum = x
			print('xnum: ',x)
		
		path = 'figures/052717/dIvsy/figVariable_'+variable+'/'
		folder = ""
		title = 'dI vs y, ' + variable + ' = ' + str(x)
		
		if variable is not 'B':
			folder+='B_%.1f_'%param.B
			title += ' B = %.1f' %param.B
		if variable is not 'Ef':
			folder+='Ef_%.1f_'%param.Ef
			title += ' Ef = %.1f' %param.Ef
		if variable is not 'ky':
			folder+='ky_%.3f_'%param.ky
			title += ' ky = %.3f' %param.ky
		if variable is not 'phi':
			folder+='phi_%.1f_'%param.phi			
			title += ' phi = %.1f' %param.phi
		if variable is not 'Bmax':
			folder+='Bmax_%.1f_'%param.Bmax			
			title += ' Bmax = %.1f' %param.Bmax
		if variable is not 'Bmin':
			folder+='Bmin_%.1f_'%param.Bmin			
			title += ' Bmin = %.1f' %param.Bmin
		if variable is not 'kyInterp':
			folder+='kyInterp_%.1f_'%param.ky_max_interp		
			title += 'kyInterp = %.1f' %param.ky_max_interp		
		if variable is not 'anum':
			folder+='anum_%.1f_'%param.anum		
			title += 'anum = %.1f' %param.anum
		if variable is not 'xnum':
			folder+='xnum_%.1f_'%param.xnum		
			title += 'xnum = %.1f' %param.xnum
		
		folder += '/'
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		name = '%s_%.2f_N_%d_n_%d_from_%.2f_to_%.2f' % (variable,x,N,n,starty,endy)
		name = name.replace('.','-')
		if param.interp:
			name+='_interp'	
	
	
	
		start_interp_time = time.time()
		obj = myObject(param)
		end_interp_time = time.time()
		time_interp = end_interp_time - start_interp_time
	
		xVal_array = np.linspace(starty,endy, N)
		F_array = np.zeros(xVal_array.shape)
		time_array = np.zeros(xVal_array.shape)
		for i in range(N):
			param.y = xVal_array[i]
			
			print('count:',i+1,'/',N)
			start = time.time()
			F_array[i] = currentDensity(param.y,param)
			end = time.time()
			time_array[i] = end-start
		
		fig = plt.figure()
		plt.plot(xVal_array,F_array,'.')		
		plt.title(title)

		fig.savefig(path+name+'.png')
		plt.close(fig)
		
		fig = plt.figure()
		plt.plot(xVal_array,time_array)
		title = 'Time vs %s for Current density. Interpolation time %.2f' %(xVariable,time_interp)
		plt.title(title)
		fig.savefig(path+name+'_time.png')
		plt.close(fig)



		
########################################################
# Current vs Phi
def plotAndSaveCurrentvsPhi(start,end,figCount,param,N):
	print('plotAndSaveCurrentvsPhi')
	B_array = np.linspace(start,end,figCount)
	init_param = copy(param)
	for B in B_array:
		#phi_start = -np.pi
		phi_start = -np.pi
		phi_end = np.pi
		phi_array = np.linspace(phi_start, phi_end, N)
		I_array = np.zeros(phi_array.shape)
		for i in range(N):
			param = init_param
			param.B = B
			param.phi = phi_array[i]
			print('count:',i+1,'/',N)
			print('      phi = ',phi_array[i])
			start = time.time()
			if param.dbl:
				I_array[i] = totalCurrent2(copy(param))
			else:
				I_array[i] = totalCurrent(copy(param))
			end = time.time()
			print('time spent: ',end-start)
			print(' ')
		fig = plt.figure()
		plt.plot(phi_array,I_array,'.')
		title = 'B = %.2f' % B
		plt.title(title)
		path = 'figures/052717/IvsPhi/'
		
		folder = ""
		folder+='Ef_%.1f_'%param.Ef
		folder+='Bmax_%.1f_'%param.Bmax			
		folder+='Bmin_%.1f_'%param.Bmin			
		folder+='kyInterp_%.1f_'%param.ky_max_interp		
		folder+='kyInt_%.1f_'%param.ky_max_int	
		folder+='anum_%.1f_'%param.anum		
		folder+='xnum_%.1f_'%param.xnum		
		folder = folder.replace('.','-')
		
		path += folder
	
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
			
		name = 'B_%.2f' % B
		name += 'W_%.f' % param.W
		name += 'y_%.f' % param.y
		if param.interp:
			name += '_interp'
		if param.dbl:
			name += '_dbl'
		name = name.replace('.','-')
	
		fig.savefig(path+name+'.png')
		plt.close(fig)

##########################################################
def plotAndSavedFreeEnergyvsy(variable,start,end,figCount,starty,endy,param,N):
#def plotAndSaveFvsPhi(start,end,figCount,B,ky,y,Ef,L,Z,N,method,variable):
	print('plotAndSaveCurrentDensityvsy')
	if variable is not 'ky':
		print('ky',param.ky)
	if variable is not 'B':
		print('B',param.B)
	if variable is not 'Ef':
		print('Ef',param.Ef)
	if variable is not 'phi':
		print('phi',param.phi)
	print('N',N)
	x_array = np.linspace(start,end,figCount)
	initial_param = copy(param)
	for x in x_array:
		param = copy(initial_param)
		if variable is 'B':
			param.B = x
			print('B: ',x)
		elif variable is 'ky':
			param.ky = x
			print('ky: ',x)
		elif variable is 'Ef':
			param.Ef = x
			print('Ef: ',x)
		elif variable is 'y':
			param.y = x
			print('y: ',x)
		elif variable is 'phi':
			param.phi = x
			print('phi: ',x)
		elif variable is 'Bmax':
			param.Bmax = x
			print('Bmax: ',x)
		elif variable is 'Bmin':
			param.Bmin = x
			print('Bmin: ',x)
		elif variable is 'kyInterp':
			param.kyInterp = x
			print('kyInterp: ',x)
		elif variable is 'anum':
			param.anum = x
			print('anum: ',x)
		elif variable is 'xnum':
			param.xnum = x
			print('xnum: ',x)
		
		path = 'figures/052717/dFvsy/figVariable_'+variable+'/'
		folder = ""
		title = 'dF vs y, ' + variable + ' = ' + str(x)
		
		if variable is not 'B':
			folder+='B_%.1f_'%param.B
			title += ' B = %.1f' %param.B
		if variable is not 'Ef':
			folder+='Ef_%.1f_'%param.Ef
			title += ' Ef = %.1f' %param.Ef
		if variable is not 'ky':
			folder+='ky_%.3f_'%param.ky
			title += ' ky = %.3f' %param.ky
		if variable is not 'phi':
			folder+='phi_%.1f_'%param.phi			
			title += ' phi = %.1f' %param.phi
		if variable is not 'Bmax':
			folder+='Bmax_%.1f_'%param.Bmax			
			title += ' Bmax = %.1f' %param.Bmax
		if variable is not 'Bmin':
			folder+='Bmin_%.1f_'%param.Bmin			
			title += ' Bmin = %.1f' %param.Bmin
		if variable is not 'kyInterp':
			folder+='kyInterp_%.1f_'%param.ky_max_interp		
			title += 'kyInterp = %.1f' %param.ky_max_interp		
		if variable is not 'anum':
			folder+='anum_%.1f_'%param.anum		
			title += 'anum = %.1f' %param.anum
		if variable is not 'xnum':
			folder+='xnum_%.1f_'%param.xnum		
			title += 'xnum = %.1f' %param.xnum
		
		
		folder += '/'
		folder = folder.replace('.','-')
		path += folder
		directory = os.path.dirname(path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		name = '%s_%.2f_N_%d_n_%d_from_%.2f_to_%.2f' % (variable,x,N,n,starty,endy)
		name = name.replace('.','-')
		if param.interp:
			name+='_interp'
		
	
		start_interp_time = time.time()
		obj = myObject(param)
		end_interp_time = time.time()
		time_interp = end_interp_time - start_interp_time
	
		xVal_array = np.linspace(starty,endy, N)
		F_array = np.zeros(xVal_array.shape)
		time_array = np.zeros(xVal_array.shape)
		for i in range(N):
			param.y = xVal_array[i]
			
			print('count:',i+1,'/',N)
			start = time.time()
			F_array[i] = dFreeEnergy(copy(param).ky,copy(param))
			end = time.time()
			time_array[i] = end-start
		
		fig = plt.figure()
		plt.plot(xVal_array,F_array,'.')
		plt.title(title)
		fig.savefig(path+name+'.png')
		plt.close(fig)
		
		fig = plt.figure()
		plt.plot(xVal_array,time_array)
		title = 'Time vs %s for dF. Interpolation time %.2f' %(xVariable,time_interp)
		plt.title(title)
		fig.savefig(path+name+'_time.png')
		plt.close(fig)

def plotAndSaveCurrentvsB(B_start,B_end,param,N):
	print('plotAndSaveCurrentvsB')
	print('k_max',param.ky_max_int)
	init_param = copy(param)
	B_array = np.linspace(B_start,B_end,N)
	phi_array = np.linspace(0.,np.pi, 10)
	I_array = np.zeros(B_array.shape)
	temp_I_array = np.zeros(phi_array.shape)
	for i in range(N):
		print('count:',i+1,'/',N)
		print('B: ',B_array[i])
		for j in range(10):
			param = copy(init_param)
			param.B = B_array[i]
			param.phi = phi_array[j]
			start = time.time()
			temp_I_array[j] = totalCurrent(param)
			end = time.time()
			print('time spent: ',end-start)
		print(' ')
		I_array[i] = np.amax(temp_I_array)
	fig = plt.figure()
	plt.plot(B_array,I_array,'.')
	
	title = 'I vs B'
	plt.title(title)
	path = 'figures/052717/IvsB/'
	
	folder = ""
	folder+='Ef_%.1f_'%param.Ef
	folder+='Bmax_%.1f_'%param.Bmax			
	folder+='Bmin_%.1f_'%param.Bmin			
	folder+='kyInterp_%.1f_'%param.ky_max_interp		
	folder+='kyInt_%.1f_'%param.ky_max_int	
	folder+='anum_%.1f_'%param.anum		
	folder+='xnum_%.1f_'%param.xnum		
	folder = folder.replace('.','-')
	
	path += folder

	directory = os.path.dirname(path)
	if not os.path.exists(directory):
		os.makedirs(directory)
		
	name = 'W_%.f' % param.W
	name += 'y_%.f' % param.y
	if param.interp:
		name += '_interp'
	name = name.replace('.','-')

	fig.savefig(path+name+'.png')
	plt.close(fig)

	
	
	
	
	
	
	title = 'Current vs magnetic field'
	plt.title(title)
	path = 'figures/052717/IvsB/'
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
	
def testFunction(B, Ef, ky, y, L, Z, N, n, method):
	print('testFunction')
	print('method',method)
	print('ky',ky)
	print('Ef',Ef)
	print('B',B)
	print('n',n)
	print('N',N)
	phi_start = -np.pi
	phi_end = np.pi
	phi = np.linspace(phi_start,phi_end,n)
	E_array = np.linspace(-1,1,N)
	Func_array = np.zeros(E_array.shape)
	plt.figure()
	for j in range(n):
		print('count:',j+1,'/',n)
		for i in range(N):		
			Func_array[i] = fun([E_array[i]],phi[j],y,B,ky,Ef,L,Z,method)[0]
		plotLabel='phi = '+str(phi[j])
		plt.plot(E_array,Func_array,label=plotLabel)
	plt.legend()
	plt.show()
	
	
def shellPlot(param):
	param.ky = 0.1
	param.B = 1.
	param.Bmin = 0
	param.Bmax = 0
	param.ky_max_interp = ky
	param.anum = 30
	param.xnum = 300
	
	N = 100
	xVariable = 'y'
	startXval = -param.L/4
	endXval = param.L/4
	
	
	B = [0.1, 1., 2., 4., 6., 8.]
	
	for i in range(0,6):
		figCount = 10
		
		variable = 'phi'
		start = -np.pi/2
		end = np.pi/2
		
		param.B = B[i]
		
		param.interp = True
		plotAndSaveFvsPhi(variable,start,end,figCount,xVariable,startXval,endXval,param,N)
		param.interp = False
		plotAndSaveFvsPhi(variable,start,end,figCount,xVariable,startXval,endXval,param,N)
	

		

ky = 0.
phi = 0.
B = 1.
Bmin = 0
Bmax = 0
ky_max = 0.
ky_max_interp = ky
anum = 10
xnum = 100
Ef = 500.
L = 106.7
W = L/4
Z = 0.
kBT = 1.
n = 4
interp = True
dbl = False
y = 0.

figCount = 10
variable = 'phi'
start = -np.pi
end = np.pi

N = 100
xVariable = 'y'
startXval = -L/4
endXval = L/4

param = parameters(y,ky,phi,B,Bmin,Bmax,ky_max,ky_max_interp,anum,xnum,Ef,L,W,Z,kBT,interp,dbl)

shellPlot(param)
#plotAndSaveEvsPhi(variable,start,end,figCount,xVariable,startXval,endXval,param,N,n)
#plotAndSaveFvsPhi(variable,start,end,figCount,xVariable,startXval,endXval,param,N)
#plotAndSaveCurrentDensityvsy(variable,start,end,figCount,-L/4,L/4,param,N)
#plotAndSavedFreeEnergyvsy(variable,start,end,figCount,-L/2,L/2,param,N)
#makePlotEvsKy(B,phi,Ef,L,Z,100,n,method)
#plotCurrentvsPhi(B,Ef,L,Z,kBT,N,method)
#testFunction(B, Ef, ky, y, L, Z, N, n, method)
#plotFvsPhi(B,Ef,ky,L,Z,kBT,N,method)
#plotAndSaveCurrentvsPhi(1.,5.,5,param,N)
#plotAndSaveCurrentvsB(0.1,5.,param,20)
#plotAndSaveCurrentvsPhi(start,end,figCount,k_max,B,Ef,L,W,Z,kBT,N,method,variable,intLim)
#plotAndSaveCurrentvsB(start,end,k_max,Ef,L,W,Z,kBT,N,method,intLim)


#####################
#To remember
#In plot from papers: Ef = 10, B = 0.5
#In produced plot: Ef=500, B = 1.0 corresponds to nu=1000 and delta=1
   



