import numpy as np

Nx = 3
Ny = 5
x = np.linspace(0,10,Nx)
y = np.linspace(30,40,Ny)
Y,X = np.meshgrid(y,x)
Z = np.zeros(X.shape)
for i in range(Nx):
	for j in range(Ny):
		Z[i,j] = X[i,j] + Y[i,j]
		
print(Z)