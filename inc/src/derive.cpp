varC integrator::derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	varC d;

	d.E = XT.E*l->g*l->sqrtkap*exp(-1.0*1.0i*l->T*l->dO)*exp(XT.G*(-0.5*1.0i*l->ag + 0.5) - XT.Q*(-0.5*1.0i*l->aq + 0.5)) - l->g*(X.E.real() + 1.0*1.0i*X.E.imag());

	d.G = -X.G*l->gg + l->Jg - (exp(X.G) - 1)*exp(-X.Q)*norm(X.E);

	d.Q = -l->rs*(-1 + exp(-X.Q))*exp(-X.Q)*norm(X.E) + (X.J + l->gq)*(-X.Q + l->q0);

	d.J = -X.J*f->wLP + f->K*f->wLP*norm(Xtau.E);


	return d;
}

