#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <omp.h>

#define M0		0.08
#define M1		0.5
#define M2		1.0
#define M3		120.0
#define G1		1.3
#define G2		2.2
#define G3		2.7
#define	RANGE	(M3 - M0)
#define ITS		10000
#define SZCL	1000000
#define STEPS	1

//#define WEIGHT(a,b,g)	(pow(a*b,-g)*(a*pow(b,g)-b*pow(a,g))/(g-1))

double func(double mass)
{
	if (mass < M1)
		return pow(0.5,-G2)/pow(0.5,-G1)*pow(mass,-G1);
	else if (mass < M2)
		return pow(mass,-G2);
	else
		return pow(mass,-G3);
}

double mass(gsl_rng *r)
{
	double Ix, pos, y, fx;
	do
	{
		Ix = 17.2601*gsl_rng_uniform_pos(r);
		pos = 0.859044*pow((17.2628-Ix),-5./6.);
		fx = func(pos);
		y = pow(pos,-G2);
		if ( fx > y ) printf("ERROR!\n");
		y *= gsl_rng_uniform_pos(r);
	} while (y > fx);
	return pos;
}

void double_print(double *arr, FILE *f)
{
	for (size_t i = 0; i < ITS; i++)
	{
		for (size_t t = 0; t < STEPS+1; t++)
			fprintf(f,"%g\t",arr[i*(STEPS+1)+t]);
		fprintf(f,"\n");
	}
	return;
}
void int_print(size_t *arr, FILE *f)
{
	for (size_t i = 0; i < ITS; i++)
	{
		for (size_t t = 0; t < STEPS+1; t++)
			fprintf(f,"%zu\t",arr[i*(STEPS+1)+t]);
		fprintf(f,"\n");
	}
	return;
}

int main(void)
{
	size_t steps = STEPS+1;
	double *sm = malloc(sizeof(double)*ITS*steps);
	double *cm = malloc(sizeof(double)*ITS*steps);
	double *lm = malloc(sizeof(double)*ITS*steps);
	size_t *ms = malloc(sizeof(size_t)*ITS*steps);
	size_t *wd = malloc(sizeof(size_t)*ITS*steps);
	size_t *ns = malloc(sizeof(size_t)*ITS*steps);
	size_t *bh = malloc(sizeof(size_t)*ITS*steps);

	#pragma omp parallel for schedule(dynamic) num_threads(16)
	for (int i = 0 ; i < ITS; i++)
	{
		// PREPARE RNG
		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_rng_env_setup();
		T = gsl_rng_default;	// Defaults to Mersenne Twister
		r = gsl_rng_alloc(T);
		gsl_rng_set(r,i);

		// construct cluster
		double *cluster = malloc(SZCL*sizeof(double));
		double *tau = malloc(SZCL*sizeof(double));
		for (size_t j = 0; j < SZCL; j++)
		{
			cluster[j] = mass(r);
			tau[j] = 1.0e10*pow(cluster[j],-2.5);
		}	

		// run time
		for (int t = 0; t < steps; t++)
		{
			size_t pos = i*steps+t;
			sm[pos] = 0;
			cm[pos] = 0;
			ms[pos] = 0;
			wd[pos] = 0;
			ns[pos] = 0;
			bh[pos] = 0;
			lm[pos] = 0;
			double lum = 0;
			double time = 10000*pow(1200000,((double) t) / STEPS);
			for (size_t s = 0; s < SZCL; s++)
			{
				if (time < tau[s])
				{
					ms[pos]++;
					sm[pos] += cluster[s];
					cm[pos] += cluster[s];
					lum += cluster[s]/tau[s];
				}
				else
				{
					if (cluster[s] < 8) // white dwarf
					{
						wd[pos]++;
						cm[pos] += 0.6;
					}
					else if (cluster[s] < 30) // neutron star
					{
						ns[pos]++;
						cm[pos] += 1.4;
					}
					else // black hole
					{
						bh[pos]++;
						cm[pos] += 10.;
					}
				}
			}
			lm[pos] = 1e-10*cm[pos] / lum;
		}

		// free memory
		gsl_rng_free(r);
		free(cluster);
		free(tau);
		printf("DONE WITH %d\n",i);
	}

	printf("PRINTING!\n");

	// print
	FILE *f;
	f = fopen("time12,log","w");
	for (int t = 0; t < steps; t++)
		fprintf(f,"%g\n",10000*pow(1500000,((double) t) / STEPS));
	fclose(f);
	f = fopen("ms12,log","w");
	int_print(ms,f);
	fclose(f);
	f = fopen("wd12,log","w");
	int_print(wd,f);
	fclose(f);
	f = fopen("ns12,log","w");
	int_print(ns,f);
	fclose(f);
	f = fopen("bh12,log","w");
	int_print(bh,f);
	fclose(f);
	f = fopen("sm12,log","w");
	double_print(sm,f);
	fclose(f);
	f = fopen("cm12,log","w");
	double_print(cm,f);
	fclose(f);
	f = fopen("lm12,log","w");
	double_print(lm,f);
	fclose(f);

	free(ms);
	free(bh);
	free(wd);
	free(lm);
	free(sm);
	free(cm);
	free(ns);

	return 0;
}
