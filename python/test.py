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


myzeros = np.zeros((3,2))
for i in range(3):
	myzeros[i][0] = i
	
print(myzeros)