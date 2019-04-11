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
#define ITS		1000
#define LIMWD	8
#define LIMNS	30
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

int exist_supernova(double thresh,int size, gsl_rng *r)
{
	int flag = 0;
	for (int i = 0; i < size; i++)
		if ( mass(r) > thresh ){flag = 1; break;}
	return flag;
}

int supernova(int size, gsl_rng *r)
{
	int sum = 0;
	for (int i = 0; i < size; i++)
		if ( mass(r) > LIMWD ) sum++;
	return sum;
}

void get_cluster(double *cluster, size_t size, gsl_rng *r)
{
	for (size_t i = 0 ; i < size; i++)
		cluster[i] = mass(r);
	return;
}

void get_lifetime(double *cluster, double *tau, size_t size)
{
	for (size_t i = 0; i < size; i++)
		tau[i] = 10*pow(cluster[i],-2.5); // in Gigayears
	return;
}

double integrate(double low, double high, size_t its, gsl_rng *r)
{
	double sum = 0;
	for (size_t i = 0; i < its; i++)
		sum += func(low+(high-low)*gsl_rng_uniform_pos(r));
	return (high-low)*((double) sum) / its;
}

int main(void)
{

	// BINS
	int binmax = 5000;
	int binmin = 50;
	int nbin = binmax - binmin;
	int counts[nbin];

	// COLLECT ALL THE CLUSTERS
	printf("Collecting clusters\n");
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0 ; i < nbin; i++)
	{
		printf("COLLECTING N = %d\n",binmin+i);
		// PREPARE RNG
		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_rng_env_setup();
		T = gsl_rng_default;	// Defaults to Mersenne Twister
		r = gsl_rng_alloc(T);
		fflush(stdout);
		counts[i] = 0;
		for (size_t j = 0 ; j < ITS; j++)
			if ( exist_supernova(8,binmin + i,r) )
				counts[i]++;
		gsl_rng_free(r);
	}

	// FIRST A
	FILE *f = fopen("supernova.dat","w");
	fprintf(f,"#SIZE\tNO. OF SUPERNOVAE\n");
	for (int i = 0 ; i < nbin; i++)
		fprintf(f,"%d\t%g\n",binmin+i,((double) counts[i])/ITS);
	fclose(f);

	// FIRST B
	int m[ITS];
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < ITS; i++)
	{
		printf("COLLECTING ITER NO. %d\n",i);
		// PREPARE RNG
		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_rng_env_setup();
		T = gsl_rng_default;	// Defaults to Mersenne Twister
		r = gsl_rng_alloc(T);
		gsl_rng_set(r,i);
		m[i] = supernova(5000,r);
		gsl_rng_free(r);
	}
	{
		FILE *f = fopen("quartiles.dat","w");
		for (int i = 0; i < ITS; i++)
		{
			fprintf(f,"%d\n",m[i]);
		}
		fclose(f);
	}

	// FIRST C
	{
		// BINS
		int binmax = 10000;
		int binmin = 50;
		int nbin = binmax - binmin;
		int counts[nbin];

		// COLLECT ALL THE CLUSTERS
		printf("Collecting clusters\n");
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0 ; i < nbin; i++)
		{
			printf("COLLECTING N = %d\n",binmin+i);
			// PREPARE RNG
			const gsl_rng_type *T;
			gsl_rng *r;
			gsl_rng_env_setup();
			T = gsl_rng_default;	// Defaults to Mersenne Twister
			r = gsl_rng_alloc(T);
			fflush(stdout);
			counts[i] = 0;
			for (size_t j = 0 ; j < ITS; j++)
				if ( exist_supernova(25,binmin + i,r) )
					counts[i]++;
			gsl_rng_free(r);
		}

		FILE *f = fopen("sun.dat","w");
		fprintf(f,"#SIZE\tNO. OF SUPERNOVAE\n");
		for (int i = 0 ; i < nbin; i++)
			fprintf(f,"%d\t%g\n",binmin+i,((double) counts[i])/ITS);
		fclose(f);
	}
	return 0;
}
