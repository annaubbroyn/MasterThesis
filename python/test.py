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

nu = 2000.

a = -nu/2
k = 0.8
L = 106.7
xi = np.sqrt(2/nu)*(L+k*nu)
res = mpmath.hyp1f1(0.5*a+0.25,0.5,0.5*xi**2)
ans = float(res.real)+1j*float(res.imag)
#ans = np.exp(-0.25*xi**2)#*
print(float(np.exp(-0.25*xi**2)*res.real))