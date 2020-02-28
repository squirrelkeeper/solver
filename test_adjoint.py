import numpy as np
import matplotlib.pyplot as pl

# dx/dt = -a * x + b * x(t-tau)

class p:
	a1 = 1.0
	a2 = 1.0
	b = 0.5

def derive(X, now, delay, p):
	dx[0] = X[now][1]
	dx[1] = -p.a1 * X[now][0]
	
	# + p.b * X[delay]
	return dx

ic_x = 1.0
dt = 1e-3
t_end = 400.0

it = int(t_end/dt)

tau = 1.0
hist = int(tau/dt)

now = 0
delay = hist

X = np.full(hist+1, ic_x)

RES = []
t = 0

for i in range(0, it):
	
	xnew = X[now] + dt * derive(X, now, delay, p)
	
	RES.append([t, xnew])
	
	
	X[delay] = xnew
	
	
	now = delay
	delay = (delay + 1)%hist
	
	
	t += dt

RES = list(zip(*RES))

pl.plot(RES[0], RES[1])
pl.show()
