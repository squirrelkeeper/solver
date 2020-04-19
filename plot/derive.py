import sympy as sy

def R(G, Q, ):
	return 0

sy.init_printing(use_unicode=True)

E = sy.Symbol('E', complex=True)
E_T = sy.Symbol('E_T', complex=True)
E_tau = sy.Symbol('E_tau', complex=True)
dE = sy.Symbol('dE', complex=True)

ER, EI, G, Q, J = sy.symbols('ER EI G Q J', real=True)
ER_T, EI_T, G_T, Q_T, J_T = sy.symbols('ER_T EI_T G_T Q_T J_T', real=True)
ER_tau, EI_tau, G_tau, Q_tau, J_tau = sy.symbols('ER_tau EI_tau G_tau Q_tau J_tau', real=True)


dER, dEI, dG, dQ, dJ = sy.symbols('dER dEI dG dQ dJ', real=True)

g, a, dO, T = sy.symbols('g a dO T', real=True)
k, ag, aq = sy.symbols('k ag aq', real=True)
Jg, gg = sy.symbols('Jg gg', real=True)
gq, q0, rs = sy.symbols('gq q0 rs', real=True)
wLP, K, tau = sy.symbols('wLP K tau', real=True)

dE = -(g+a*1j)*E
dER = sy.re(dE)
print(dER)





dER = -g * ER + a * EI
dEI = -g * EI - a * ER

dG = Jg - gg * G - sy.exp(-Q)*(sy.exp(G) - 1)*(ER*ER+EI*EI)
dQ = (gq + J) * (q0 - Q) - rs * sy.exp(-Q) * (sy.exp(-Q) - 1)*(ER*ER+EI*EI)
dJ = -wLP * J

dER+= 4







Psi = sy.Matrix([ER, EI, G, Q, J])
dPsi = sy.Matrix([dER, dEI, dG, dQ, dJ])

A = dPsi.jacobian(Psi)
APsi = A*Psi
