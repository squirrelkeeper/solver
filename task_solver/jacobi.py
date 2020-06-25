import sympy as sy
import numpy as np
import matplotlib.pyplot as pl


g, gg, gq = sy.symbols('g gg gq')
ag, aq = sy.symbols('ag aq')
Jg, q0 =  sy.symbols('Jg q0')
rs, sqrtkap = sy.symbols('rs sqrtkap')


K, tau, wLP = sy.symbols('K tau wLP')

dX = sy.Matrix(sy.symbols('dX0:5'))
X = sy.Matrix(sy.symbols('X0:5'))
XT = sy.Matrix(sy.symbols('XT0:5'))
Xtau = sy.Matrix(sy.symbols('Xtau0:5'))

dY = sy.Matrix(sy.symbols('dY0:5')).T
Y = sy.Matrix(sy.symbols('Y0:5')).T
YT = sy.Matrix(sy.symbols('YT0:5')).T
Ytau = sy.Matrix(sy.symbols('Ytau0:5')).T

Z = sy.Matrix(sy.symbols('Z0:5'))

dX[0] =-g * X[0]
dX[0]+= g * sqrtkap * sy.exp((XT[2]-XT[3])/2) * sy.cos((aq * XT[3] - ag * XT[2])/2) * XT[0]
dX[0]-= g * sqrtkap * sy.exp((XT[2]-XT[3])/2) * sy.sin((aq * XT[3] - ag * XT[2])/2) * XT[1]

dX[1] =-g * X[1]
dX[1]+= g * sqrtkap * sy.exp((XT[2]-XT[3])/2) * sy.cos((aq * XT[3] - ag * XT[2])/2) * XT[1]
dX[1]+= g * sqrtkap * sy.exp((XT[2]-XT[3])/2) * sy.sin((aq * XT[3] - ag * XT[2])/2) * XT[0]

dX[2] = Jg - gg * X[2] - sy.exp(-X[3]) * (sy.exp(X[2])-1) * (X[0]**2+X[1]**2)

dX[3] = (gq + X[4])*(q0 - X[3]) - rs * sy.exp(-X[3]) * (sy.exp(X[3])-1) * (X[0]**2+X[1]**2)

dX[4] = wLP * K * (Xtau[0]**2+Xtau[1]**2) - wLP * X[4]


A = dX.jacobian(X)
B = dX.jacobian(XT)
C = dX.jacobian(Xtau)


dY = Y*A + YT*B + Ytau*C

b11, b12, b13, b14, b21, b22, b23, b24 = sy.symbols('b11 b12 b13 b14 b21 b22 b23 b24')
B = sy.Matrix([[b11, b12, b13, b14, 0],[b21, b22, b23, b24, 0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])


c51, c52 = sy.symbols('c51 c52')
C = sy.Matrix([[0, 0, 0, 0, 0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[c51,c52,0,0,0]])

Z[2,0] = 0
Z[3,0] = 0
Z[4,0] = 0
#print(YT*B*Z)
print(Y*C*Z)
