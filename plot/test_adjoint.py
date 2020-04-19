import numpy as np
import matplotlib.pyplot as pl

# dx/dt = -a * x + b * x(t-tau)

class p:
	a = 0.0
	K = 0.2
	D = 0.0
	J = 3.0
	T = 1.0

def derive(X, tp, p, t):
	dx = np.zeros(3)
	
	dx[0] = X[tp[0]][0] * X[tp[0]][2]
	dx[0]+= -p.a * X[tp[0]][1] * X[tp[0]][2]
	dx[0]+= p.K * X[tp[1]][0]
	dx[0]+= p.D
	
	dx[1] = X[tp[0]][1] * X[tp[0]][2]
	dx[1]+=  p.a * X[tp[0]][0] * X[tp[0]][2]
	dx[1]+= p.K * X[tp[1]][1]
	dx[1]+= p.D
	
	dx[2] = (p.J - X[tp[0]][2] - (X[tp[0]][1] + 1.0) * (X[tp[0]][0]*X[tp[0]][0]+X[tp[0]][1]*X[tp[0]][1]))/p.T
	
	return dx

def derive_inearized(X, Y, tp, p, t):









def hist_init(X, ic_X, opt):
	if opt == "const":
		for i in range(0, len(X)):
			for j in range(0, len(X[0])):
				X[i][j] = ic_X[j]

ic_x0 = 1.0
ic_x1 = 1.0
ic_x2 = 1.0
ic_X = [ic_x0, ic_x1, ic_x2]

dt = 1e-3
t_end = 400.0
t_out = 400.0


it = int(t_end/dt)

tau = 1.0
hist = int(tau/dt)

t = 0.0
tp = [0, hist]

X = np.zeros((hist+1, 3))
hist_init(X, ic_X, "const")


RES = np.zeros((int(t_out/dt), 4))

for i in range(0, it):
	
	xnew = X[tp[0]] + dt * derive(X, tp, p, t)
	
	X[tp[1]] = xnew
	
	tp[0] = tp[1]
	tp[1] = (tp[1] + 1)%hist
	
	t += dt

	xnew = np.insert(xnew, 0, t)
	

	RES[i] = xnew




RES = np.asarray(list(zip(*RES)))


#pl.plot(RES[1], RES[2])

I = RES[1]*RES[1]+RES[2]*RES[2]


pl.plot(RES[0], I)
pl.plot(RES[0], RES[3])

pl.show()
