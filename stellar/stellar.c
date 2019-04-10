#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define M0		0.08
#define M1		0.5
#define M2		1.0
#define M3		120.0
#define G1		1.3
#define G2		2.2
#define G3		2.7
#define	RANGE	(M3 - M0)
#define THRESH	8
#define ITS		1000000

double mass(gsl_rng *r)
{
	double top = pow(M0,-G1) > pow(M1,-G2) ? pow(M0,-G1) : pow(M1,-G2);
	double fx, y;
	do
	{
		double pos = M0 + RANGE*gsl_rng_uniform_pos(r);
		if (pos > M2) fx = pow(pos,-G3);
		else if (pos > M1) fx = pow(pos,-G2);
		else fx = pow(pos,-G1);
		y = top*gsl_rng_uniform_pos(r);
	} while (y > fx)
	return fx;
}

int supernova(int size, gsl_rng *r)
{
	int sum = 0;
	for (int i = 0; i < size; i++)
		sum += mass(r) > THRESH ? 1 : 0;
	return sum;
}



int main(void)
{
	binmax = 5000;
	binmin = 50;
	nbin = binmax - binmin;
	size_t counts[nbin];

	for (size_t i = 0 ; i < nbin; i++)
	{
		counts[i] = 0;
		for (size_t j = 0 ; j < ITS; j++)
			if ( supernova(binmin + i,r) )
				counts[i]++;
	}

	// FIRST A
	FILE *f = fopen(supernova.dat);
	fprintf("#SIZE\tNO. OF SUPERNOVAE\n");
	for (size_t i = 0 ; i < nbin; i++)
		fprintf(f,"%zu\t%zu\n",binmin+i,counts[i]);

	// FIRST B
	gsl_rstat *rstat = gsl_rstat_alloc();
	for (size_t i = 0; i < ITS; i++)
		gsl_rstat_add(supernova(5000,r));
	printf("MEAN: %g\n",gsl_rstat_mean(rstat));
	printf("MEDIAN: %g\n",gsl_rstat_median(rstat));
	printf("UPPER Q.: %g\n",gsl_rstat_mean(rstat));
	printf("LOWER Q.: %g\n",gsl_rstat_mean(rstat));
	gsl_rstat_free(rstat);

	// Follow-up

	return 0;
}
