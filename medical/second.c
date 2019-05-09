#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_histogram.h>
#include <omp.h>

#define NBIN	500
#define NITS	10000000

void compton(double ine, double *sce, double *sca, gsl_rng *r)
{
	double lamb = ine/0.5109989;
	double q, rho, r1, r2, r3, check;
	int track;
	q = (2.*lamb+1.) / (2.*lamb+9.);
	do
	{
		r1 = gsl_rng_uniform_pos(r);
		r2 = gsl_rng_uniform_pos(r);
		r3 = gsl_rng_uniform_pos(r);
		if ( r1 < q )
		{
			track = 1;
			rho = 1.+2.*lamb*r2;
			check = 4.*(rho-1.)/(rho*rho);
		}
		else
		{
			track = 2;
			rho = (2.*lamb+1.)/(2.*lamb*r2+1);
			check = 1-(rho-1)/lamb;
			check = 0.5*(check*check + 1./rho);
		}
	} while (r3 > check);
	*sce = ine/rho;

	*sca = track == 1 ? acos(1-2*r2) : acos(1-(rho-1)/lamb);
}


int main(void)
{
	// PREPARE RNG and HISTOGRAM
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_histogram *hist1 = gsl_histogram_alloc(NBIN);
	gsl_histogram *hist2 = gsl_histogram_alloc(NBIN);
	gsl_histogram_set_ranges_uniform(hist1,0,M_PI);
	gsl_histogram_set_ranges_uniform(hist2,0,1);
	
	// FILES and ACCUMULATOR
	FILE *f1, *f2, *f3;
	f1 = fopen("secondi.dat","w");
	f2 = fopen("secondii.dat","w");
	int nj = 90;
	double bot[nj];
	double top[nj];
	double de = 0.9 / nj;

	for (int j = 0; j < nj; j++)
	{
		double e = 0.1 + de * j;
		gsl_histogram_reset(hist1);
		gsl_histogram_reset(hist2);
		printf("e: %g\n",e);
		for (size_t i = 0; i < NITS; i++)
		{
			double sca, sce;
			compton(e,&sce,&sca,r);
			gsl_histogram_increment(hist1,sca);
			gsl_histogram_increment(hist2,1-sce/e);
		}
		fprintf(f1,"%g",e);
		fprintf(f2,"%g",e);
		double max = 0;
		for (int i = 0; i < NBIN; i++)
		{
			fprintf(f1,"\t%g",hist1->bin[i]);
			fprintf(f2,"\t%g",hist2->bin[i]);
			if (hist2->bin[i] > max)
			{
				max = hist2->bin[i];
				bot[j] = hist2->range[i];
				top[j] = hist2->range[i+1];
			}
		}
		fprintf(f1,"\n");
		fprintf(f2,"\n");
	}
	fclose(f1);
	fclose(f2);

	gsl_histogram_free(hist1);
	gsl_histogram_free(hist2);

	f3 = fopen("pairs.dat","w");
	for (size_t j = 0; j < nj; j++)
	{
		double e = 0.1 + de * j;
		printf("e: %g\n",e);
		double angle = 0;
		size_t count = 0;
		for (size_t i = 0; i < NITS; i++)
		{
			double sca, sce;
			compton(e,&sce,&sca,r);
			if ((1-sce/e) >= bot[j] && (1-sce/e) < top[j])
			{
				count++;
				angle += sca;
			}
		}
		printf("count: %zu\n",count);
		fprintf(f3,"%g\t%g\n",e,angle/count/M_PI*180);
	}
	fclose(f3);

	gsl_rng_free(r);

	return 0;
}
