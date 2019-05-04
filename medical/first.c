#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rstat.h>
#include <omp.h>



double get_path(double mu, gsl_rng *r)
{
	return -log(gsl_rng_uniform_pos(r))/mu;
}

void paths(double mu, gls_rstat_workspace* mean, gsl_rstat_quantile_workspace* low, gsl_rstat_quantile_workspace* high, gsl_rng *r)
{
	gsl_rstat_quantile_reset(low);
	gsl_rstat_quantile_reset(high);
	gsl_rstat_reset(rstat);

	for (size_t i = 0; i < NITS; i++)
	{
		double d = get_path(mu,r);
		gsl_rstat_add(rstat);
		gsl_rstat_quantile_add(low);
		gsl_rstat_quantile_add(high);
	}

	return
}

void doit(double rho, FILE *fin, FILE *fout)
{
	gsl_rstat_workspace *rstat = gsl_rstat_alloc();
	gsl_rstat_quantile_workspace *low = gsl_rstat_quantile_alloc(0.25);
	gsl_rstat_quantile_workspace *high = gsl_rstat_quantile_alloc(0.75);

	char line[100];
	double cs;
	double e;

	rewind(fin);
	while (fgets(line,sizeof line, fin) != NULL)
	{
		sscanf(line,"%lf\t%lf",&e,&cs);
		paths(*rho,rstat,low,high,r);
		fprintf(fout,"%g\t%g\t%g\t%g\t%g\n",
			e,gsl_rstat_mean(rstat),gsl_rstat_median(rstat),
			gsl_rstat_quantile_get(low),gsl_rstat_quantile_get(high));
	}

	gsl_rstat_free(rstat);
	gsl_rstat_quantile_free(low);
	gsl_rstat_quantile_free(high);
}

int main(void)
{
	// PREPARE RNG
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	
	double rho;
	FILE *fin, *fout;

	// water
	fin = fopen("water.dat","r");
	fout = fopen("water.out","w");
	rho = 1;
	doit(rho,fin,fout);
	fclose(fin);
	fclose(fout);

	// lead
	fin = fopen("lead.dat","r");
	fout = fopen("lead.out","w");
	rho = 1;
	doit(rho,fin,fout);
	fclose(fin);
	fclose(fout);

	// aluminium
	fin = fopen("al.dat","r");
	fout = fopen("al.out","w");
	rho = 1;
	doit(rho,fin,fout);
	fclose(fin);
	fclose(fout);

	// iodine
	fin = fopen("iodine.dat","r");
	fout = fopen("iodine.out","w");
	rho = 1;
	doit(rho,fin,fout);
	fclose(fin);
	fclose(fout);

	gsl_rng_free(r);

	return 0;
}
