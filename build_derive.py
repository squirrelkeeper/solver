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

########################################################################

class par:
	g, dw, T = sy.symbols('l->g l->dw l->T', positive=True)
	sqrtkap, ag, aq = sy.symbols('l->sqrtkap l->ag l->aq', positive=True)
	Jg, gg = sy.symbols('l->Jg l->gg', positive=True)
	gq, q0, rs = sy.symbols('l->gq l->q0 l->rs', positive=True)
	wLP, K, tau = sy.symbols('f->wLP f->K f->tau', positive=True)

########################################################################

def R(G_T, Q_T,  p):
	R = p.g * p.sqrtkap
	R*= sy.exp(0.5*(1 - 1j * p.ag)*G_T - 0.5*(1 - 1j * p.aq)*Q_T)
	R*= sy.exp(-1j * p.dw * p.T)
	return R

def RR(G_T, Q_T, p):
	RR = p.g * p.sqrtkap
	RR*= sy.exp(0.5 * G_T - 0.5 * Q_T)
	RR*= sy.cos(0.5 * (p.aq * Q_T - p.ag * G_T) - p.dw * p.T)
	return RR
	
def RI(G_T, Q_T, p):
	RI = p.g * p.sqrtkap
	RI*= sy.exp(0.5 * G_T - 0.5 * Q_T)
	RI*= sy.sin(0.5 * (p.aq * Q_T - p.ag * G_T) - p.dw * p.T)
	return RI

########################################################################

def derive_full():
	p = par()
	
	dE, E, E_T, E_tau = sy.symbols('d.E X.E XT.E Xtau.E', complex=True)

	dG, G, G_T, G_tau = sy.symbols('d.G X.G XT.G Xtau.G', real=True)
	dQ, Q, Q_T, Q_tau = sy.symbols('d.Q X.Q XT.Q Xtau.Q', real=True)
	dJ, J, J_T, J_tau = sy.symbols('d.J X.J XT.J Xtau.J', real=True)
	
	dE = -p.g * (sy.re(E) + sy.I * sy.im(E)) + R(G_T, Q_T, p) * E_T
	
	dG = p.Jg - p.gg * G - sy.exp(-Q)*(sy.exp(G) - 1)*sy.Abs(E)**2
	dQ = (p.gq + J) * (p.q0 - Q) - p.rs * sy.exp(-Q) * (sy.exp(Q) - 1)*sy.Abs(E)**2
	dJ = p.wLP * (p.K * sy.Abs(E_tau)**2 - J)

	dE = sy.powsimp(sy.expand(dE))
	dE = sy.expand(dE, complex=True)
	dG = sy.powsimp(sy.expand(dG))
	dQ = sy.powsimp(sy.expand(dQ))
	dJ = sy.powsimp(sy.expand(dJ))
	
	Psi = sy.Matrix([E, G, Q, J])
	dPsi = sy.Matrix([dE, dG, dQ, dJ])
	
	coll = [Psi, dPsi]
	
	return coll

def build_derive_full(coll):
	PsiC = coll[0]
	dPsiC = coll[1]
	rule = '\nvarC integrator::derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)'
	rule+='\n'
	rule+= '{'
	rule+='\n\t'
	rule+='varC d;'
	rule+='\n\n'
	for i in range(0,len(dPsiC)):
		eq = '\t'
		eq+= 'd.'
		eq+= str(PsiC[i])[2::]
		eq+=' = '
		eq+= str(dPsiC[i])
		eq+=';\n\n'
		
		eq = eq.replace('exp(', 'expf(')
		eq = eq.replace('sin(', 'sinf(')
		eq = eq.replace('cos(', 'cosf(')
		
		eq = eq.replace('Abs(X.E)**2', 'norm(X.E)')
		eq = eq.replace('Abs(Xtau.E)**2', 'norm(Xtau.E)')

		eq = eq.replace('re(X.E)', 'X.E.real()')
		eq = eq.replace('re(XT.E)', 'XT.E.real()')
		eq = eq.replace('im(X.E)', 'X.E.imag()')
		eq = eq.replace('im(XT.E)', 'XT.E.imag()')
		
		eq = repl_I(eq)
		eq = repl_square(eq)
		
		rule+= eq
	
	rule+='\n\treturn d;'
	rule+='\n'
	rule+= '}'

	return rule

########################################################################

def derive_ret():
	p = par()
	
	dER, dEI, dG, dQ, dJ = sy.symbols('d.ER d.EI d.G d.Q d.J', real=True)
	
	ER, ER_T, ER_tau = sy.symbols('X.ER XT.ER Xtau.ER', real=True)
	EI, EI_T, EI_tau = sy.symbols('X.EI XT.EI Xtau.EI', real=True)
	G, G_T, G_tau = sy.symbols('X.G XT.G Xtau.G', real=True)
	Q, Q_T, Q_tau = sy.symbols('X.Q XT.Q Xtau.Q', real=True)
	J, J_T, J_tau = sy.symbols('X.J XT.J Xtau.J', real=True)

	ERr, ERr_T, ERr_tau = sy.symbols('Y.ER YT.ER Ytau.ER', real=True)
	EIr, EIr_T, EIr_tau = sy.symbols('Y.EI YT.EI Ytau.EI', real=True)
	Gr, Gr_T, Gr_tau = sy.symbols('Y.G YT.G Ytau.G', real=True)
	Qr, Qr_T, Qr_tau = sy.symbols('Y.Q YT.Q Ytau.Q', real=True)
	Jr, Jr_T, Jr_tau = sy.symbols('Y.J YT.J Ytau.J', real=True)

	
	dER = -p.g * ER
	dER+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.cos(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * ER_T
	dER+= -p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.sin(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * EI_T
	
	
	dEI = -p.g * EI
	dEI+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.cos(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * EI_T
	dEI+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.sin(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * ER_T
	
	dG = p.Jg - p.gg * G 
	dG-= sy.exp(-Q)*(sy.exp(G) - 1) * (ER*ER+EI*EI)
	
	dQ = (p.gq + J) * (p.q0 - Q) 
	dQ-= p.rs * sy.exp(-Q) * (sy.exp(Q) - 1) * (ER*ER+EI*EI)
	
	dJ = p.wLP * (p.K * (ER_tau*ER_tau+EI_tau*EI_tau) - J)

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
	
	for i in range(0,len(dPsir)):
		dPsir[i] = sy.powsimp(sy.expand(dPsir[i]))
	
	coll = [Psir, dPsir]
		
	return coll

def build_derive_ret(coll):
	Psir = coll[0]
	dPsir = coll[1]
	rule = '\nvar integrator::derive_ret(var &X, var &XT, var &Xtau, var &Y, var &YT, var &Ytau, lpar_dbl_set *l, fpar_dbl_set *f)'
	rule+='\n'
	rule+= '{'
	rule+= '\n\t//X is the homogenous solution, Y the pertubation'
	rule+='\n\t'
	rule+='var d;'
	rule+='\n\n'
	for i in range(0,len(dPsir)):
		eq = '\t'
		eq+= 'd.'
		eq+= str(Psir[i])[2::]
		eq+=' = '
		eq+= str(dPsir[i])
		eq+=';\n\n'
		
		eq = eq.replace('exp(', 'expf(')
		eq = eq.replace('sin(', 'sinf(')
		eq = eq.replace('cos(', 'cosf(')
		
		eq = repl_square(eq)
		
		rule+= eq
	
	rule+='\n\treturn d;'
	rule+='\n'
	rule+= '}'

	return rule

########################################################################

def derive_adj():
	p = par()
	
	dER, dEI, dG, dQ, dJ = sy.symbols('d.ER d.EI d.G d.Q d.J', real=True)
	
	ER, ER_T, ER_tau = sy.symbols('X.ER XT.ER Xtau.ER', real=True)
	EI, EI_T, EI_tau = sy.symbols('X.EI XT.EI Xtau.EI', real=True)
	G, G_T, G_tau = sy.symbols('X.G XT.G Xtau.G', real=True)
	Q, Q_T, Q_tau = sy.symbols('X.Q XT.Q Xtau.Q', real=True)
	J, J_T, J_tau = sy.symbols('X.J XT.J Xtau.J', real=True)

	ERa, ERa_T, ERa_tau = sy.symbols('Z.ER ZT.ER Ztau.ER', real=True)
	EIa, EIa_T, EIa_tau = sy.symbols('Z.EI ZT.EI Ztau.EI', real=True)
	Ga, Ga_T, Ga_tau = sy.symbols('Z.G ZT.G Ztau.G', real=True)
	Qa, Qa_T, Qa_tau = sy.symbols('Z.Q ZT.Q Ztau.Q', real=True)
	Ja, Ja_T, Ja_tau = sy.symbols('Z.J ZT.J Ztau.J', real=True)

	
	dER = -p.g * ER
	dER+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.cos(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * ER_T
	dER+= -p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.sin(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * EI_T
	
	
	dEI = -p.g * EI
	dEI+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.cos(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * EI_T
	dEI+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.sin(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * ER_T
	
	dG = p.Jg - p.gg * G 
	dG-= sy.exp(-Q)*(sy.exp(G) - 1) * (ER*ER+EI*EI)
	
	dQ = (p.gq + J) * (p.q0 - Q) 
	dQ-= p.rs * sy.exp(-Q) * (sy.exp(Q) - 1) * (ER*ER+EI*EI)
	
	dJ = p.wLP * (p.K * (ER_tau*ER_tau+EI_tau*EI_tau) - J)


	Psi = sy.Matrix([ER, EI, G, Q, J])
	Psi_T = sy.Matrix([ER_T, EI_T, G_T, Q_T, J_T])
	Psi_tau = sy.Matrix([ER_tau, EI_tau, G_tau, Q_tau, J_tau])

	Psia = sy.Matrix([ERa, EIa, Ga, Qa, Ja])
	Psia_T = sy.Matrix([ERa_T, EIa_T, Ga_T, Qa_T, Ja_T])
	Psia_tau = sy.Matrix([ERa_tau, EIa_tau, Ga_tau, Qa_tau, Ja_tau])

	dPsi = sy.Matrix([dER, dEI, dG, dQ, dJ])

	A = dPsi.jacobian(Psi)
	B = dPsi.jacobian(Psi_T)
	C = dPsi.jacobian(Psi_tau)

	dPsia = Psia.T * A + Psia_T.T * B + Psia_tau.T * C
	
	for i in range(0,len(dPsia)):
		dPsia[i] = sy.powsimp(sy.expand(dPsia[i]))
	
	coll = [Psia, dPsia]
		
	return coll

def build_derive_adj(coll):
	Psia = coll[0]
	dPsia = coll[1]
	rule = '\nvar integrator::derive_adj(var &X, var &XT, var &Xtau, var &Z, var &ZT, var &Ztau, lpar_dbl_set *l, fpar_dbl_set *f)'
	rule+='\n'
	rule+= '{'
	rule+= '\n\t//X is the homogenous solution, Z the adjoint pertubation'
	rule+='\n\t'
	rule+='var d;'
	rule+='\n\n'
	for i in range(0,len(dPsia)):
		eq = '\t'
		eq+= 'd.'
		eq+= str(Psia[i])[2::]
		eq+=' = '
		eq+= str(dPsia[i])
		eq+=';\n\n'
		
		eq = eq.replace('exp(', 'expf(')
		eq = eq.replace('sin(', 'sinf(')
		eq = eq.replace('cos(', 'cosf(')
		
		eq = repl_square(eq)
		
		rule+= eq
	
	rule+='\n\treturn d;'
	rule+='\n'
	rule+= '}'

	return rule

########################################################################

def bilinear():
	p = par()

	dER, dEI, dG, dQ, dJ = sy.symbols('d.ER d.EI d.G d.Q d.J', real=True)

	ER, ER_T, ER_tau = sy.symbols('X.ER XT.ER Xtau.ER', real=True)
	EI, EI_T, EI_tau = sy.symbols('X.EI XT.EI Xtau.EI', real=True)
	G, G_T, G_tau = sy.symbols('X.G XT.G Xtau.G', real=True)
	Q, Q_T, Q_tau = sy.symbols('X.Q XT.Q Xtau.Q', real=True)
	J, J_T, J_tau = sy.symbols('X.J XT.J Xtau.J', real=True)

	ERr, ERr_T, ERr_tau = sy.symbols('Y.ER YT.ER Ytau.ER', real=True)
	EIr, EIr_T, EIr_tau = sy.symbols('Y.EI YT.EI Ytau.EI', real=True)
	Gr, Gr_T, Gr_tau = sy.symbols('Y.G YT.G Ytau.G', real=True)
	Qr, Qr_T, Qr_tau = sy.symbols('Y.Q YT.Q Ytau.Q', real=True)
	Jr, Jr_T, Jr_tau = sy.symbols('Y.J YT.J Ytau.J', real=True)

	ERa, ERa_T, ERa_tau = sy.symbols('Z.ER ZT.ER Ztau.ER', real=True)
	EIa, EIa_T, EIa_tau = sy.symbols('Z.EI ZT.EI Ztau.EI', real=True)
	Ga, Ga_T, Ga_tau = sy.symbols('Z.G ZT.G Ztau.G', real=True)
	Qa, Qa_T, Qa_tau = sy.symbols('Z.Q ZT.Q Ztau.Q', real=True)
	Ja, Ja_T, Ja_tau = sy.symbols('Z.J ZT.J Ztau.J', real=True)


	dE = -p.g * (ER + sy.I * EI) + (RR(G_T, Q_T, p) + sy.I * RI(G_T, Q_T, p)) * (ER_T + sy.I * EI_T)

	dER = sy.re(dE)
	dEI = sy.im(dE)

	dG = p.Jg - p.gg * G - sy.exp(-Q)*(sy.exp(G) - 1)*(ER**2+EI**2)
	dQ = (p.gq + J) * (p.q0 - Q) - p.rs * sy.exp(-Q) * (sy.exp(Q) - 1)*(ER**2+EI**2)
	dJ = p.wLP * (p.K * (ER_tau**2+EI_tau**2) - J)

	Psi = sy.Matrix([ER, EI, G, Q, J])
	Psi_T = sy.Matrix([ER_T, EI_T, G_T, Q_T, J_T])
	Psi_tau = sy.Matrix([ER_tau, EI_tau, G_tau, Q_tau, J_tau])

	Psir = sy.Matrix([ERr, EIr, Gr, Qr, Jr])
	Psir_T = sy.Matrix([ERr_T, EIr_T, Gr_T, Qr_T, Jr_T])
	Psir_tau = sy.Matrix([ERr_tau, EIr_tau, Gr_tau, Qr_tau, Jr_tau])

	Psia = sy.Matrix([ERa, EIa, Ga, Qa, Ja])
	Psia_T = sy.Matrix([ERa_T, EIa_T, Ga_T, Qa_T, Ja_T])
	Psia_tau = sy.Matrix([ERa_tau, EIa_tau, Ga_tau, Qa_tau, Ja_tau])

	dPsi = sy.Matrix([dER, dEI, dG, dQ, dJ])

	B = dPsi.jacobian(Psi_T)
	C = dPsi.jacobian(Psi_tau)

	Asum = Psia.T * Psir
	Bsum = Psia_T.T * B * Psir_T
	Csum = Psia_tau.T * C * Psir_tau

	coll = [Asum,Bsum,Csum]
	
	return coll

def build_bilinear(coll):
	eq = [str(coll[0][0]),str(coll[1][0]),str(coll[2][0])]
	
	for i in range(0,len(eq)):
		eq[i] = eq[i].replace('exp(', 'expf(')
		eq[i] = eq[i].replace('sin(', 'sinf(')
		eq[i] = eq[i].replace('cos(', 'cosf(')
		
		
		eq[i] = eq[i].replace('X.', 'X[pos0].')	
		eq[i] = eq[i].replace('Y.', 'Y[pos0].')
		eq[i] = eq[i].replace('Z.', 'Z[pos0].')
		
		eq[i] = eq[i].replace('XT.', 'X[pos0+r].')		
		eq[i] = eq[i].replace('YT.', 'Y[pos0+r].')
		eq[i] = eq[i].replace('ZT.', 'Z[pos0+r+dim1].')
		
		eq[i] = eq[i].replace('Xtau.', 'X[pos0+r].')
		eq[i] = eq[i].replace('Ytau.', 'Y[pos0+r].')
		eq[i] = eq[i].replace('Ztau.', 'Z[pos0+r+dim2].')	
	
	
	A = eq[0]
	B = eq[1]
	C = eq[2]
	
	rule = '\ndouble integrator::bilinear_step(vector<var> &X, vector<var> &Y, vector<var> &Z, lpar_dbl_set *l, fpar_dbl_set *f)'
	rule+='\n'
	rule+= '{'
	rule+= '\n\t//X is the homogenous solution, Y the pertubation, Z the adjoint pertubation'
	rule+='\n\t'

	rule+='\n\n'
	rule+='\tdouble b = '+A+';\n\n'
	
	rule+='\tfor(long r = -dim1; r < 0; r++)\n'
	rule+='\t{\n'
	
	rule+='\t\t'
	rule+='b+='+B+';\n'
	rule+='\t}'
	rule+='\n\n'
	
	
	rule+='\tfor(long r = -dim2; r < 0; r++)\n'
	rule+='\t{\n'
	rule+='\t\t'
	rule+='b+='+C+';\n'
	rule+='\t}\n'
	
	rule+='\n\treturn b;'
	rule+='\n'
	rule+= '}'

	return rule

########################################################################

def derive_real():
	p = par()
	
	dER, ER, ER_T, ER_tau = sy.symbols('d.ER X.ER XT.ER Xtau.ER', real=True)
	dEI, EI, EI_T, EI_tau = sy.symbols('d.EI X.EI XT.EI Xtau.EI', real=True)

	dG, G, G_T, G_tau = sy.symbols('d.G X.G XT.G Xtau.G', real=True)
	dQ, Q, Q_T, Q_tau = sy.symbols('d.Q X.Q XT.Q Xtau.Q', real=True)
	dJ, J, J_T, J_tau = sy.symbols('d.J X.J XT.J Xtau.J', real=True)
	
	dER = -p.g * ER
	dER+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.cos(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * ER_T
	dER+= -p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.sin(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * EI_T
	
	
	dEI = -p.g * EI
	dEI+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.cos(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * EI_T
	dEI+=  p.g * p.sqrtkap * sy.exp(0.5*(G_T-Q_T)) * sy.sin(0.5*(p.aq * Q_T - p.ag * G_T) - p.T * p.dw) * ER_T
	
	dG = p.Jg - p.gg * G 
	dG-= sy.exp(-Q)*(sy.exp(G) - 1) * (ER*ER+EI*EI)
	
	dQ = (p.gq + J) * (p.q0 - Q) 
	dQ-= p.rs * sy.exp(-Q) * (sy.exp(Q) - 1) * (ER*ER+EI*EI)
	
	dJ = p.wLP * (p.K * (ER_tau*ER_tau+EI_tau*EI_tau) - J)

	dER = sy.powsimp(sy.expand(dER))
	dEI = sy.powsimp(sy.expand(dEI))
	dG = sy.powsimp(sy.expand(dG))
	dQ = sy.powsimp(sy.expand(dQ))
	dJ = sy.powsimp(sy.expand(dJ))
	
	Psi = sy.Matrix([ER, EI, G, Q, J])
	dPsi = sy.Matrix([dER, dEI, dG, dQ, dJ])
	
	coll = [Psi, dPsi]
	
	return coll

def build_derive_real(coll):
	PsiC = coll[0]
	dPsiC = coll[1]
	rule = '\n'
	rule+='var integrator::derive_real(var &X, var &XT, var &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)'
	rule+='\n'
	rule+= '{'
	rule+='\n\t'
	rule+='var d;'
	rule+='\n\n'
	for i in range(0,len(dPsiC)):
		eq = '\t'
		eq+= 'd.'
		eq+= str(PsiC[i])[2::]
		eq+=' = '
		eq+= str(dPsiC[i])
		eq+=';\n\n'
		
		eq = eq.replace('exp(', 'expf(')
		eq = eq.replace('sin(', 'sinf(')
		eq = eq.replace('cos(', 'cosf(')

		eq = repl_square(eq)
		
		rule+= eq
	
	rule+='\n\treturn d;'
	rule+='\n'
	rule+= '}'

	return rule









sy.init_printing(use_unicode=True)











all_func = ''


all_func+= "/*REPLACE START*/"

all_func+= build_derive_real(derive_real())
all_func+= '\n\n'
#all_func+= build_derive_full(derive_full())
#all_func+= '\n'
all_func+= build_derive_ret(derive_ret())
all_func+= '\n\n'
all_func+= build_derive_adj(derive_adj())
all_func+= '\n\n'
all_func+= build_bilinear(bilinear())
all_func+= "/*REPLACE END*/"






file_cpp = open("inc/src/integrate.cpp", "rt")
file_txt = file_cpp.read()

repl_begin = file_txt.find("/*REPLACE START*/")
repl_end = file_txt.find("/*REPLACE END*/")+18

#print(file_txt[repl_begin:repl_end:])

#to_replace = "/*REPLACE START*/"
to_replace = file_txt[repl_begin:repl_end:]
#to_replace+= "/*REPLACE END*/"

file_txt = file_txt.replace(to_replace, all_func)
file_cpp.close()

file_cpp = open("inc/src/integrate.cpp", "wt")
file_cpp.write(file_txt)
file_cpp.close()

