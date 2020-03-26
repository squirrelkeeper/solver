import sympy as sy

class par:
	g, a, dO, T = sy.symbols('g a dO T', positive=True)
	k, ag, aq = sy.symbols('k ag aq', positive=True)
	Jg, gg = sy.symbols('Jg gg', positive=True)
	gq, q0, rs = sy.symbols('gq q0 rs', positive=True)
	wLP, K, tau = sy.symbols('wLP K tau', positive=True)

def R(G_T, Q_T,  p):
	R = p.g * sy.sqrt(p.k)
	R*= sy.exp(0.5*(1 - 1j * p.ag)*G_T - 0.5*(1 - 1j * p.aq)*Q_T)
	R*= sy.exp(-1j*(p.dO + p.a) * p.T)
	return R

def RR(G_T, Q_T, p):
	RR = p.g * sy.sqrt(p.k)
	RR*= sy.exp(0.5 * G_T - 0.5 * Q_T)
	RR*= sy.cos(0.5 * (p.aq * Q_T - p.ag * G_T) - (p.dO + p.a) * p.T)
	return RR
	
def RI(G_T, Q_T, p):
	RI = p.g * sy.sqrt(p.k)
	RI*= sy.exp(0.5 * G_T - 0.5 * Q_T)
	RI*= sy.sin(0.5 * (p.aq * Q_T - p.ag * G_T) - (p.dO + p.a) * p.T)
	return RI

sy.init_printing(use_unicode=True)

p = par()

dE, E, E_T, E_tau = sy.symbols('dE E E_T E_tau', complex=True)

dER, ER, ER_T, ER_tau = sy.symbols('dER ER ER_T ER_tau', real=True)
dEI, EI, EI_T, EI_tau = sy.symbols('dEI EI EI_T EI_tau', real=True)

dG, G, G_T, G_tau = sy.symbols('dG G G_T G_tau', real=True)
dQ, Q, Q_T, Q_tau = sy.symbols('dQ Q Q_T Q_tau', real=True)
dJ, J, J_T, J_tau = sy.symbols('dJ J J_T J_tau', real=True)


dE = -(p.g + p.a * 1j) * (ER + 1j * EI)

dER = sy.re(dE)
dER+= RR(G_T, Q_T, p) * ER_T - RI(G_T, Q_T, p) * EI_T

dEI = sy.im(dE)
dEI+= RR(G_T, Q_T, p) * EI_T + RI(G_T, Q_T, p) * ER_T


dG = p.Jg - p.gg * G - sy.exp(-Q)*(sy.exp(G) - 1)*(ER*ER+EI*EI)
dQ = (p.gq + J) * (p.q0 - Q) - p.rs * sy.exp(-Q) * (sy.exp(-Q) - 1)*(ER*ER+EI*EI)
dJ = -p.wLP * J
dJ+= p.wLP * p.K * (ER_tau*ER_tau+EI_tau*EI_tau)


Psi = sy.Matrix([ER, EI, G, Q, J])
Psi_T = sy.Matrix([ER_T, EI_T, G_T, Q_T, J_T])
Psi_tau = sy.Matrix([ER_tau, EI_tau, G_tau, Q_tau, J_tau])

dPsi = sy.Matrix([dER, dEI, dG, dQ, dJ])

A = dPsi.jacobian(Psi)
B = dPsi.jacobian(Psi_T)
C = dPsi.jacobian(Psi_tau)

dPsi_ret = A * Psi + B * Psi_T + C * Psi_tau

dPsi_adj = Psi.T * A + Psi_T.T * B + Psi_tau.T * C

BIL_1 = Psi.T * Psi
BIL_2 = Psi.T
BIL_3 = 

print(dPsi_adj[0])




