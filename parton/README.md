# Parton Showers

*I have doubts about the obtained event numbers and thus have either a bug in the code or a misunderstanding in my interpretation of how the veto algorithm works.*

We consider the Parton shower with the processes
`q->qg` and `g->gg`.
These have the rate function
```
fabc(p,z) = alpha_s / (2*pi) Pabc(z) / p^2
```
for `p < E_a / 2` and zero otherwise.

The spatial distributions Pabc(z) are
`Pqqg(z) = 4/3*(1+z^2)/(1-z)` and `Pggg(z)=3(1-z(1-z))^2/(z*(1-z))`.

We note that both `Pqqg` and `Pggg` are maximized in the limit `z=1-1/E_a`
and use veto sampling with
```
gabc(p,z) = alpha_s / (2*pi) Pabc(1-1/E_a) / p^2 / (1-2 p/E_a)
```
Integrated over the `z`-coordinate from `p/E_a` to `1-p/E_a` we have
```
gabc(p) = alpha_s / (2*pi) Pabc(1-1/E_a) / p^2
```
The veto algorithm is implemented in the file `parton.c` as
```C
void veto(double ea, double pmax, double (*pabc)(double z), gsl_rng* r, double *pout, double *zout)
{
	double term, z, y, g, ft, p, p2;

	y = 1;
	ft = 0;
	p2 = pmax*pmax;
	term = 0;
	do {
		// determine time
		term = term - log(gsl_rng_uniform_pos(r));
		p2 = pmax*pmax*exp(-2*M_PI*term/(ALPHA*pabc(1-1/ea)));
		p = sqrt(p2);
		if ( 2*p < ea )
		{
			// determine space
			z = p/ea + gsl_rng_uniform_pos(r)*(1.-2.*p/ea);
			g = 0.5*M_1_PI*ALPHA*pabc(1-1/ea)/(p2*(1-2*p/ea));
			y = gsl_rng_uniform_pos(r)*g;
			ft = 0.5*M_1_PI*ALPHA/p2*pabc(z);
		}
	} while (y > ft);
	*pout = p;
	*zout = z;

	return;
}
```
The shower is handled by calling the veto algorithm recursively for the
two daughter particles with new energies and new `pmax=p`,
and counting the number of `q->qg` events as well as the number of added gluons.
This is implemented in the file `parton.c` as the functions `runq` for the
quark and `rung` for the gluon.
```C
void runq(double ea, double pmax, size_t *s, gsl_rng *r)
{
	if (ea <= 1 || pmax <= 1) return;

	// do veto and find outcome
	double p, z, eb, ec;
	veto(ea,pmax,&q2qg,r,&p,&z);
	eb = ea*z;
	ec = ea-eb;

	s[0]++; // q->qg event
	s[1]++; // gluon created

	// recursion
	runq(eb,p,s,r);
	rung(ec,p,s,r);

	return;
}
```
with similar procedure for the gluon.

Sampling 100 thousand times we obtain

|`E_0`		| `avg(q->qg)`	| `sd(q->qg)`	| `avg(#g)`	| `sd(#g)`|
|-----------|---------------|---------------|-----------|---------|
|50|1.80367|0.821006|3.30272|2.65958|
|100|2.16755|0.987991|5.038|4.02471|
|200|2.6073|1.15288|7.93786|6.1498|
|400|3.10205|1.30999|12.7074|9.43829|

The above numbers seem low. But at least the means increase with increasing energy.
I had expected a Poissonian process and hence Poissonian SDs, which I have not obtained.
As I cannot recognise any bugs in the code, I must have misinterpreted some step in the algorithm.
