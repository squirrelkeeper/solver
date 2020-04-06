import sympy as sy
import string

def var(s):
	if s.isalnum() or s=='_' or s=='.':
		return True
	else:
		return False

def repl_square(expr):
	ind = []
	r1 = []
	r2 =[]
	for i in range(0,len(expr)):
		if expr[i]=='*' and expr[i-1]=='*' and i > 1:
			ind.append(i-2)
	for i in range(0, len(ind)):
		b = ''
		p = 0
		for j in range(ind[i],0,-1):
			if var(expr[j]) and p < 1:
				b+= expr[j]
			elif expr[j]==')' and p < 1 and j == ind[i]:
				p = 1
				b+= expr[j]
			elif p == 1 and not expr[j]=='(':
				b+= expr[j]
			elif p == 1 and expr[j]=='(':
				b+= expr[j]
				p = 0
			else:
				break
		b = b[::-1]
		r2.append(b)
		t = '#'+ str(i) +'#'
		expr = expr[0:ind[i]+1:] + t + expr[ind[i]+4::]
		r1.append(t)	
	for i in range(0,len(r1)):
		expr = expr.replace(r1[i],'*'+r2[i])
	return expr


def repl_I(expr):
	expr = list(expr)
	
	for i in range(0, len(expr)):
		if expr[i]=='I' and not var(expr[i-1]) and not var(expr[i+1]):
			expr[i] = '1.0i'
	
	expr = ''.join(expr)
	return expr

class par:
	g, dO, T = sy.symbols('l->g l->dO l->T', positive=True)
	sqrtkap, ag, aq = sy.symbols('l->sqrtkap l->ag l->aq', positive=True)
	Jg, gg = sy.symbols('l->Jg l->gg', positive=True)
	gq, q0, rs = sy.symbols('l->gq l->q0 l->rs', positive=True)
	wLP, K, tau = sy.symbols('f->wLP f->K f->tau', positive=True)

def R(G_T, Q_T,  p):
	R = p.g * p.sqrtkap
	R*= sy.exp(0.5*(1 - 1j * p.ag)*G_T - 0.5*(1 - 1j * p.aq)*Q_T)
	R*= sy.exp(-1j * p.dO * p.T)
	return R

def RR(G_T, Q_T, p):
	RR = p.g * p.sqrtkap
	RR*= sy.exp(0.5 * G_T - 0.5 * Q_T)
	RR*= sy.cos(0.5 * (p.aq * Q_T - p.ag * G_T) - p.dO * p.T)
	return RR
	
def RI(G_T, Q_T, p):
	RI = p.g * p.sqrtkap
	RI*= sy.exp(0.5 * G_T - 0.5 * Q_T)
	RI*= sy.sin(0.5 * (p.aq * Q_T - p.ag * G_T) - p.dO * p.T)
	return RI

def build_derive_full(coll):
	PsiC = coll[0]
	dPsiC = coll[1]
	derive_rule_full = 'varC integrator::derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)'
	derive_rule_full+='\n'
	derive_rule_full+= '{'
	derive_rule_full+='\n\t'
	derive_rule_full+='varC d;'
	derive_rule_full+='\n\n'
	for i in range(0,len(dPsiC)):
		full_eq = '\t'
		full_eq+= 'd.'
		full_eq+= str(PsiC[i])[2::]
		full_eq+=' = '
		full_eq+= str(dPsiC[i])
		full_eq+=';\n\n'
		
		full_eq = full_eq.replace('Abs(X.E)**2', 'norm(X.E)')
		full_eq = full_eq.replace('Abs(Xtau.E)**2', 'norm(Xtau.E)')

		full_eq = full_eq.replace('re(X.E)', 'X.E.real()')
		full_eq = full_eq.replace('im(X.E)', 'X.E.imag()')
		
		full_eq = repl_I(full_eq)
		full_eq = repl_square(full_eq)
		
		derive_rule_full+= full_eq
	
	derive_rule_full+='\n\treturn d;'
	derive_rule_full+='\n'
	derive_rule_full+= '}'

	return derive_rule_full

def derive_full():
	p = par()
	
	dE, E, E_T, E_tau = sy.symbols('d.E X.E XT.E Xtau.E', complex=True)

	dG, G, G_T, G_tau = sy.symbols('d.G X.G XT.G Xtau.G', real=True)
	dQ, Q, Q_T, Q_tau = sy.symbols('d.Q X.Q XT.Q Xtau.Q', real=True)
	dJ, J, J_T, J_tau = sy.symbols('d.J X.J XT.J Xtau.J', real=True)
	
	dE = -p.g * (sy.re(E) + 1.0j * sy.im(E))
	dE+= R(G_T, Q_T, p) * E_T
	
	dG = p.Jg - p.gg * G - sy.exp(-Q)*(sy.exp(G) - 1)*sy.Abs(E)**2
	dQ = (p.gq + J) * (p.q0 - Q) - p.rs * sy.exp(-Q) * (sy.exp(-Q) - 1)*sy.Abs(E)**2
	dJ = -p.wLP * J
	dJ+= p.wLP * p.K * sy.Abs(E_tau)**2
	
	
	Psi = sy.Matrix([E, G, Q, J])
	dPsi = sy.Matrix([dE, dG, dQ, dJ])
	
	coll = [Psi, dPsi]
	
	return coll
	

sy.init_printing(use_unicode=True)

p = par()

dER, ER, ER_T, ER_tau= sy.symbols('d.ER X.ER XT.ER Xtau.ER', real=True)
ERr, ERr_T, ERr_tau= sy.symbols('Y.ER YT.ER Ytau.ER', real=True)

dEI, EI, EI_T, EI_tau= sy.symbols('d.EI X.EI XT.EI Xtau.EI', real=True)
EIr, EIr_T, EIr_tau= sy.symbols('Y.EI YT.EI Ytau.EI', real=True)

dG, G, G_T, G_tau= sy.symbols('d.G X.G XT.G Xtau.G', real=True)
Gr, Gr_T, Gr_tau= sy.symbols('Y.G YT.G Ytau.G', real=True)

dQ, Q, Q_T, Q_tau= sy.symbols('d.Q X.Q XT.Q Xtau.Q', real=True)
Qr, Qr_T, Qr_tau= sy.symbols('Y.Q YT.Q Ytau.Q', real=True)

dJ, J, J_T, J_tau= sy.symbols('d.J X.J XT.J Xtau.J', real=True)
Jr, Jr_T, Jr_tau= sy.symbols('Y.J YT.J Ytau.J', real=True)


dER = sy.re(-p.g * (ER + 1.0j * EI))
dER+= RR(G_T, Q_T, p) * ER_T - RI(G_T, Q_T, p) * EI_T

dEI = sy.im(-p.g * (ER + 1.0j * EI))
dEI+= RR(G_T, Q_T, p) * EI_T + RI(G_T, Q_T, p) * ER_T

dG = p.Jg - p.gg * G - sy.exp(-Q)*(sy.exp(G) - 1)*(ER*ER+EI*EI)
dQ = (p.gq + J) * (p.q0 - Q) - p.rs * sy.exp(-Q) * (sy.exp(-Q) - 1)*(ER*ER+EI*EI)
dJ = -p.wLP * J
dJ+= p.wLP * p.K * (ER_tau*ER_tau+EI_tau*EI_tau)


Psi = sy.Matrix([ER, EI, G, Q, J])
Psi_T = sy.Matrix([ER_T, EI_T, G_T, Q_T, J_T])
Psi_tau = sy.Matrix([ER_tau, EI_tau, G_tau, Q_tau, J_tau])

Psir = sy.Matrix([ERr, EIr, Gr, Qr, Jr])
Psir_T = sy.Matrix([ERr_T, EIr_T, Gr_T, Qr_T, Jr_T])
Psir_tau = sy.Matrix([ERr_tau, EIr_tau, Gr_tau, Qr_tau, Jr_tau])

dPsi = sy.Matrix([dER, dEI, dG, dQ, dJ])

A = dPsi.jacobian(Psi)
B = dPsi.jacobian(Psi_T)
C = dPsi.jacobian(Psi_tau)

dPsir = A * Psir + B * Psir_T + C * Psir_tau


'''


A = dPsi.jacobian(Psi)
B = dPsi.jacobian(Psi_T)
C = dPsi.jacobian(Psi_tau)

dPsir = A * Psir + B * Psir_T + C * Psir_tau

dPsia = Psi.T * A + Psi_T.T * B + Psi_tau.T * C

BIL_1 = Psi.T * Psi
BIL_2 = Psi.T
BIL_3 = 6


'''




derive_rule_ret = ''
'''
derive_rule_full = 'varC integrator::derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)'
derive_rule_full+='\n'
derive_rule_full+= '{'
derive_rule_full+='\n\t'
derive_rule_full+='varC dX;'
derive_rule_full+='\n\n'

for i in range(0,len(dPsi)):
	full_eq = '\t'
	full_eq+= 'dX.'
	full_eq+= str(Psir[i])
	full_eq+=' = '
	full_eq+= str(dPsir[i])
	full_eq+=';\n\n'
	
#	full_eq = full_eq.replace('Abs(E)**2', 'norm(E)')
#	full_eq = full_eq.replace('Abs(E_tau)**2', 'norm(E_tau)')

#	full_eq = full_eq.replace('re(E)', 'E.real()')
#	full_eq = full_eq.replace('im(E)', 'E.imag()')
	
#	full_eq = repl_I(full_eq)
#	full_eq = repl_square(full_eq)
	
	derive_rule_ret+= full_eq

derive_rule_full+='\n\treturn dX;'
derive_rule_full+='\n'
derive_rule_full+= '}'
'''









func_derive_full = build_derive_full(derive_full())


#print(func_derive_full)

file_cpp = open("derive.cpp", "w")
file_cpp.write(func_derive_full)
file_cpp.close()



