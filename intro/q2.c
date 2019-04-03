#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

void mci_cosflat(gsl_rng* r, gsl_rstat_workspace *rstat)
{
	double R, fx;
	// select x in (0,1)
	R = gsl_rng_uniform_pos(r);	// in interval (0,1)
	// collect integral
	fx = cos(0.5*M_PI*R);
	gsl_rstat_add(fx,rstat);	// dx = dGinv(x)
	return;
}

void mci_cosimp(gsl_rng* r, gsl_rstat_workspace *rstat)
{
	double R, x, fx;
	// select select x acc. to g(x) = 1-x^2
	R = gsl_rng_uniform_pos(r);
	x = 2*cos(acos(-R)/3-2*M_PI/3);
	// collect integral
	fx = cos(0.5*M_PI*x);
	gsl_rstat_add(2./3.*fx/(1-x*x),rstat);	// dx = dGinv(x)
	return;
}

int main(void)
{
	// prepare RNG
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;	// Defaults to Mersenne Twister
	r = gsl_rng_alloc(T);

	// prepare running statistics
	gsl_rstat_workspace *flat = gsl_rstat_alloc();
	gsl_rstat_workspace *imp = gsl_rstat_alloc();

	// MC integration
	size_t tot = 1000000;			// desired no. of points	
	FILE *f = fopen("q2.out","w");	// output file
	for (size_t i = 0; i < tot; i++)
	{
		mci_cosflat(r,flat);
		mci_cosflat(r,imp);
		fprintf(f,"%zu\t%g\t%g\n",
			i,gsl_rstat_mean(flat),gsl_rstat_mean(imp));
	}
	fclose(f);
	
	// free memory
	gsl_rng_free(r);
	gsl_rstat_free(flat);
	gsl_rstat_free(imp);

	return 0;
}
