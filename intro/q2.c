#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>

// hit and miss integration, flat g(x) = 1, 0 < x < 1
double hmi_cosflat(gsl_rng* r)
{
	double R, y, fx;
	// select x in (0,1)
	R = gsl_rng_uniform_pos(r);	// in interval (0,1)
	// collect integral
	y = gsl_rng_uniform_pos(r);
	fx = cos(0.5*M_PI*R);
	return y < fx ? 1 : 0;
}

// importance sampling with g(x) = 1 - x^2
double hmi_cosimp(gsl_rng* r)
{
	double R, x, y, fx;
	// select select x acc. to g(x)
	R = gsl_rng_uniform_pos(r);
	x = 2*cos(acos(-R)/3-2*M_PI/3);
	// collect integral
	y = gsl_rng_uniform_pos(r)*(1-x*x);
	fx = cos(0.5*M_PI*x);
	return y < fx ? 2./3. : 0 ; // G(1) = 2/3
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
		gsl_rstat_add(hmi_cosflat(r),flat);
		gsl_rstat_add(hmi_cosimp(r),imp);
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
