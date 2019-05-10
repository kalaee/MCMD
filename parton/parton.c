#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rstat.h>
#include "/nfs/users2/kalaee/Documents/core/core_rstat.h"

#define REPS	1000000
#define ALPHA	0.2
#define GF(t,func)		(0.5*M_1_PI*ALPHA*func(1-1/ea)/(t))
#define GINT(t,func)	(0.5*M_1_PI*ALPHA*func(1-1/ea)*log((pmax*pmax)/(t)))
#define GINV(g,func)	((pmax*pmax)*exp(-2*M_PI*(g)/(ALPHA*func(1-1/ea))))

// http://home.thep.lu.se/~leif/MCMD/lec3.pdf

// spatial function for q -> qg
double q2qg(double z)
{
	return 4./3.*(1+z*z)/(1-z);
}

// spatial function for g -> gg
double g2gg(double z)
{
	double term = z*(1-z);
	return 3*(1-term)*(1-term)/term;
}

void veto(double ea, double pmax, double (*pabc)(double z), gsl_rng* r, double *pout, double *zout)
{
	double term, z, y, ft, p, p2;

	// determine time
	y = 1;
	ft = 0;
	p2 = pmax*pmax;
	do {
		term = GINT(p2,pabc)-log(gsl_rng_uniform_pos(r));
		p2 = GINV(term,pabc);
		p = sqrt(p2);
		if (2*p < ea)
		{
			z = p/ea + gsl_rng_uniform_pos(r)*(1.-2.*p/ea);
			y = gsl_rng_uniform_pos(r)*GF(p2,pabc);
			ft = 0.5*M_1_PI*ALPHA/p2*pabc(z);
		}
	} while (y > ft);
	*pout = p;
	*zout = z;

	return;
}

// time evolution for gluon
void rung(double ea, double pmax, size_t *s, gsl_rng *r)
{
	if (ea <= 1 || pmax <= 1) return;
	// do veto and find outcome
	double p, z, eb, ec;
	veto(ea,pmax,&g2gg,r,&p,&z);
	eb = ea*z;
	ec = ea-eb;

	s[1]++; // gluon created

	// recursion
	rung(eb,p,s,r);
	rung(ec,p,s,r);

	return;
}

// time evolution for quark
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

int main(void)
{
	// PREPARE RNG
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// other preparations
	double ea[] = {50, 100, 200, 400};
	int len = sizeof(ea) / sizeof(ea[0]);
	FILE *f = fopen("parton.log","w");
	size_t s[2];
	core_rstat *events = core_rstat_alloc(); // running statistics
	core_rstat *gluons = core_rstat_alloc(); // running statistics

	// fun over different initial configs
	for (int i = 0; i < len; i++)
	{
		core_rstat_reset(events);
		core_rstat_reset(gluons);
		for (size_t j = 0; j < REPS; j++) // samples showers
		{
			s[0] = 0;
			s[1] = 0;
			runq(ea[i],ea[i],s,r);
			core_rstat_add(s[0],events);
			core_rstat_add(s[1],gluons);
		}
		printf("%g\t%g\t%g\t%g\t%g\n",ea[i],
			core_rstat_mean(events),core_rstat_sd(events),
			core_rstat_mean(gluons),core_rstat_sd(gluons));
		fprintf(f,"%g\t%g\t%g\t%g\t%g\n",ea[i],
			core_rstat_mean(events),core_rstat_sd(events),
			core_rstat_mean(gluons),core_rstat_sd(gluons));
	}
	fclose(f);

	// free memory
	core_rstat_free(events);
	core_rstat_free(gluons);
	gsl_rng_free(r);

	return 0;
}
