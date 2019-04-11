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
#define THRESH	8
#define ITS		10000
#define LIMWD	8
#define LIMNS	30
#define WEIGHT(a,b,g)	(pow(a*b,-g)*(a*pow(b,g)-b*pow(a,g))/(g-1))

double mass(gsl_rng *r)
{
	double mass[] = {M0, M1, M2, M3};
	double gamm[] = {G1, G2, G3};
	double w0 = WEIGHT(mass[0],mass[1],gamm[0]);
	double w1 = WEIGHT(mass[1],mass[2],gamm[1]);
	double w2 = WEIGHT(mass[2],mass[3],gamm[2]);
	double sector = (w0+w1+w2)*gsl_rng_uniform_pos(r);
	int flag;
	if (sector < w0) flag = 0;
	else if (sector < w1) flag = 1;
	else flag = 2;
	double top, pos, y, fx;
	do
	{
		top = pow(mass[flag],-gamm[flag]);
		pos = mass[flag]+(mass[flag+1]-mass[flag])*gsl_rng_uniform_pos(r);
		y = top*gsl_rng_uniform_pos(r);
		fx = pow(pos,-gamm[flag]);
	} while (y > fx);
	return pos;
}

int exist_supernova(int size, gsl_rng *r)
{
	int flag = 0;
	for (int i = 0; i < size; i++)
		if ( mass(r) > THRESH ){flag = 1; break;}
	return flag;
}

int supernova(int size, gsl_rng *r)
{
	int sum = 0;
	for (int i = 0; i < size; i++)
		if ( mass(r) > THRESH ) sum++;
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

int main(void)
{

	// BINS
	int binmax = 5000;
	int binmin = 50;
	int nbin = binmax - binmin;
	int counts[nbin];

//	// COLLECT ALL THE CLUSTERS
//	printf("Collecting clusters\n");
//	#pragma omp parallel for schedule(dynamic) num_threads(20)
//	for (int i = 0 ; i < nbin; i++)
//	{
//		printf("COLLECTING N = %d\n",i);
//		// PREPARE RNG
//		const gsl_rng_type *T;
//		gsl_rng *r;
//		gsl_rng_env_setup();
//		T = gsl_rng_default;	// Defaults to Mersenne Twister
//		r = gsl_rng_alloc(T);
//		fflush(stdout);
//		counts[i] = 0;
//		for (size_t j = 0 ; j < ITS; j++)
//			if ( exist_supernova(binmin + i,r) )
//				counts[i]++;
//		gsl_rng_free(r);
//	}
//
//	// FIRST A
//	FILE *f = fopen("supernova.dat","w");
//	fprintf(f,"#SIZE\tNO. OF SUPERNOVAE\n");
//	for (int i = 0 ; i < nbin; i++)
//		fprintf(f,"%d\t%g\n",binmin+i,((double) counts[i])/ITS);
//	fclose(f);

	// FIRST B
	gsl_rstat_workspace *rstat = gsl_rstat_alloc();
	gsl_rstat_quantile_workspace *qlow = gsl_rstat_quantile_alloc(0.25);
	gsl_rstat_quantile_workspace *qhigh = gsl_rstat_quantile_alloc(0.75);
	double m[ITS];
	#pragma omp parallel for schedule(dynamic) num_threads(20)
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

		double m[i] = supernova(5000,r);

		gsl_rng_free(r);
	}
	for (int i = 0; i < ITS; i++)
	{
		gsl_rstat_add(m[i],rstat);
		gsl_rstat_quantile_add(m[i],qlow);
		gsl_rstat_quantile_add(m[i],qhigh);
	}
	FILE *f = fopen("quartiles.out","w");
	fprintf(f,"MEAN: %g\n",gsl_rstat_mean(rstat));
	fprintf(f,"MEDIAN: %g\n",gsl_rstat_median(rstat));
	fprintf(f,"LOWER Q.: %g\n",gsl_rstat_quantile_get(qlow));
	fprintf(f,"UPPER Q.: %g\n",gsl_rstat_quantile_get(qhigh));
	fclose(f);
	gsl_rstat_free(rstat);
	gsl_rstat_quantile_free(qlow);
	gsl_rstat_quantile_free(qhigh);

	// FOLLOW-UP


	return 0;
}
