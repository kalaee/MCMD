#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
//#include <gsl/gsl_rstat.h>
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
#define LIMWD	8
#define LIMNS	30
#define SZCL	1000000
#define STEPS	15

//#define WEIGHT(a,b,g)	(pow(a*b,-g)*(a*pow(b,g)-b*pow(a,g))/(g-1))

double func(double mass)
{
	if (mass < M1)
		return pow(mass,-G1);
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
		pos = 0.859044/pow((17.2628-Ix),5./6.);
		fx = func(pos);
		y = gsl_rng_uniform_pos(r)*pow(pos,-G2);
	} while (y > fx);
	return pos;
}

void double_print(double arr[ITS][STEPS+1], FILE *f)
{
	for (int i = 0; i < ITS; i++)
	{
		for (int t = 0; t < STEPS+1; t++)
			fprintf(f,"%g\t",arr[i][t]);
		fprintf(f,"\n");
	}
	return;
}
void int_print(size_t arr[ITS][STEPS+1], FILE *f)
{
	for (int i = 0; i < ITS; i++)
	{
		for (int t = 0; t < STEPS+1; t++)
			fprintf(f,"%zu\t",arr[i][t]);
		fprintf(f,"\n");
	}
	return;
}

int main(void)
{
	int steps = STEPS+1;
	double sm[ITS][steps];	// stellar mass
	double lm[ITS][steps];	// mass/luminosity ratio
	size_t ms[ITS][steps];	// main sequence
	size_t wd[ITS][steps];	// white dwarfs
	size_t ns[ITS][steps];	// neutron star
	size_t bh[ITS][steps];	// black hole

	#pragma omp parallel for schedule(dynamic) num_threads(22)
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
			tau[j] = 10*pow(cluster[j],-2.5);
		}	

		// run time
		for (int t = 0; t < steps; t++)
		{
			sm[i][t] = 0;
			ms[i][t] = 0;
			wd[i][t] = 0;
			ns[i][t] = 0;
			bh[i][t] = 0;
			lm[i][t] = 0;
			double lum = 0;
			for (size_t s = 0; s < SZCL; s++)
			{
				if (t < tau[s])
				{
					ms[i][t]++;
					sm[i][t] += cluster[s];
					lum += tau[s]*cluster[s];
				}
				else
				{
					if (cluster[s] < LIMWD) // white dwarf
					{
						wd[i][t]++;
						sm[i][t] += 0.6;
					}
					else if (cluster[s] < LIMNS) // neutron star
					{
						ns[i][t]++;
						sm[i][t] += 1.4;
					}
					else // black hole
					{
						bh[i][t]++;
						sm[i][t] += 10.;
					}
				}
			}
			lm[i][t] = lum / sm[i][t];
		}

		// free memory
		gsl_rng_free(r);
		free(cluster);
		free(tau);
		printf("DONE WITH %d\n",i);
	}

	// print
	FILE *f;
	f = fopen("ms.log","w");
	int_print(ms,f);
	fclose(f);
	f = fopen("wd.log","w");
	int_print(wd,f);
	fclose(f);
	f = fopen("ns.log","w");
	int_print(ns,f);
	fclose(f);
	f = fopen("bh.log","w");
	int_print(bh,f);
	fclose(f);
	f = fopen("sm.log","w");
	double_print(sm,f);
	fclose(f);
	f = fopen("lm.log","w");
	double_print(lm,f);
	fclose(f);

	return 0;
}
