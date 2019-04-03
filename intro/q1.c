#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

double rng_sinc(gsl_rng* r, gsl_rstat_workspace *rstat)
{
	double R, y, x, fx;
	do
	{
		// select x acc. to g(x) = 1 (x < 1) + 1/x^2 (x >= 1)
		R = 2*gsl_rng_uniform_pos(r);	// in interval (0,1)-->(0,2)
		x = R < 1 ? R : 1/(2-R);		// invert x = Ginv(R)
		// compare y = R*g(x) and f(x)
		y = gsl_rng_uniform_pos(r);
		if ( x > 1 ) y = y/(x*x);
		fx = gsl_pow_2(sin(x)/x);
		// collect integral
		if ( x > 1 ) gsl_rstat_add(fx*x*x,rstat);	// dx = x^2 * dGinv(x)
		else gsl_rstat_add(fx,rstat);				// dx = dGinv(x)
	} while (y > fx);
	return x;
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
	gsl_rstat_workspace *rstat = gsl_rstat_alloc();
	gsl_rstat_reset(rstat);

	size_t tot = 10000000;			// desired no. of points	
	FILE *f = fopen("q1.out","w");	// output file
	for (size_t i = 0; i < tot; i++)
		fprintf(f,"%g\n",rng_sinc(r,rstat));
	fclose(f);
	double volume = 2;				// Integration volume, u in (0,2)
	printf("INTEGRAL: %g\n",volume*gsl_rstat_mean(rstat));
	size_t n = gsl_rstat_n(rstat);
	printf("UNCERTAINTY: %g\n",volume*gsl_rstat_sd(rstat)/sqrt(n));
	printf("FCT CALLS: %zu\n",n);
	
	// free memory
	gsl_rng_free(r);
	gsl_rstat_free(rstat);

	return 0;
}
