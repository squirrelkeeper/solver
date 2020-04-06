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



test_string = 'dX.Q = -l->rs*(-1 + exp(-Q))*exp(-Q)*Abs(E)**2 + (J + l->gq)*(-Q + l->q0) * (1+2)**2;'
print(test_string)


test_string = repl_square(test_string)

test_string = repl_I(test_string)

print(test_string)

