#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rstat.h>
#include <omp.h>

#define NITS	1000000

double get_path(double mu, gsl_rng *r)
{
	return -log(gsl_rng_uniform_pos(r))/mu;
}

void paths(double mu, gsl_rstat_workspace* mean, gsl_rstat_quantile_workspace* low, gsl_rstat_quantile_workspace* high, gsl_rng *r)
{
	for (size_t i = 0; i < NITS; i++)
	{
		double d = get_path(mu,r);
		gsl_rstat_add(d,mean);
		gsl_rstat_quantile_add(d,low);
		gsl_rstat_quantile_add(d,high);
	}

	return;
}

void doit(double rho, FILE *fin, FILE *fout, gsl_rng *r)
{
	char line[100];
	double cs;
	double e;

	rewind(fin);
	while (fgets(line,sizeof line, fin) != NULL)
	{
		sscanf(line,"%lf\t%lf",&e,&cs);
		printf("ENERGY: %g\n",e);
		gsl_rstat_workspace *rstat = gsl_rstat_alloc();
		gsl_rstat_quantile_workspace *low = gsl_rstat_quantile_alloc(0.25);
		gsl_rstat_quantile_workspace *high = gsl_rstat_quantile_alloc(0.75);
		paths(cs*rho,rstat,low,high,r);
		fprintf(fout,"%g\t%g\t%g\t%g\t%g\n",
			e,gsl_rstat_mean(rstat),gsl_rstat_median(rstat),
			gsl_rstat_quantile_get(low),gsl_rstat_quantile_get(high));
		gsl_rstat_free(rstat);
		gsl_rstat_quantile_free(low);
		gsl_rstat_quantile_free(high);
	}

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
	doit(rho,fin,fout,r);
	fclose(fin);
	fclose(fout);

	// lead
	fin = fopen("lead.dat","r");
	fout = fopen("lead.out","w");
	rho = 11.34;
	doit(rho,fin,fout,r);
	fclose(fin);
	fclose(fout);

	// aluminium
	fin = fopen("al.dat","r");
	fout = fopen("al.out","w");
	rho = 2.7;
	doit(rho,fin,fout,r);
	fclose(fin);
	fclose(fout);

	// iodine
	fin = fopen("iodine.dat","r");
	fout = fopen("iodine.out","w");
	rho = 4.94;
	doit(rho,fin,fout,r);
	fclose(fin);
	fclose(fout);

	gsl_rng_free(r);

	return 0;
}
