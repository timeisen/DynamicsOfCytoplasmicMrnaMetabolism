/* kinetics_simple.c__________________________________________________________*/
#include <R.h>
#include <Rmath.h>

static double p[4*1+4]; /* the number of parameters, genes*3 right not */
/* Function to initialize the parameters:.....................................*/
void ode_p_init(void (* odep)(int *, double *))
{
	int N=4*1+4	; /* number of parameters again */
	odep(&N, p);
}
/* Function to compute derivatives of the ODE:................................*/
void ode_deriv (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
	if (ip[0] <1) error("nout should be at least 1");
	/* Assign static definitions: */
	#define G   	   1   						/* Number of genes   		  */
	#define b		   p[4*G+0]							/* Decapping scale param      */
	#define size       p[4*G+1]							/* tx size param for nbinom   */
	#define location   p[4*G+2]
	#define scale      p[4*G+3]
	#define T   	   251							/* max number of tails   */
	#define cutoff     101
	/* Pre-assign the rate constant arrays:							  */
	int g;
	int l;
	double k_a;
	double st_v; 
	double tx;   
	double k_dea_1; /* k_dea_1 acts first on longer tails */
	double k_dea_2;
	double k_dea;
	double k_dec;

	/* iteration and derivative assignment							  */


	for ( g = 0; g < G; g++) /* genes first	  */
	{
		st_v  = p[4*g+0]; /* starting tail length */
		k_a   = p[4*g+1]; /* transcription rate */
		k_dea_1 = p[4*g+2]; /* deadenylation rate */
		k_dea_2 = p[4*g+3]; /* deadenylation rate */
		k_dea = k_dea_2+plogis(0,cutoff,1,1,0)*(k_dea_1-k_dea_2);
		tx    = dnbinom_mu(T, size, st_v, 0)*k_a;
		k_dec = plogis(0,location,scale,1,0)*b; /* decay rate */
		ydot[g*T] = tx - (k_dea + k_dec)*y[g*T];
		for (l = 1; l < (T-1); l++)
		{
			tx   = dnbinom_mu(T - l, size, st_v, 0)*k_a;
			k_dea = k_dea_1+plogis(l,cutoff,1,1,0)*(k_dea_2-k_dea_1);
			k_dec  = plogis(l,location,scale,1,0)*b;
			ydot[l + g*T] = tx - (k_dea + k_dec)*y[l + g*T] + k_dea*y[l - 1 + g*T];
		}
		k_dea = k_dea_1+plogis(T-1,cutoff,1,1,0)*(k_dea_2-k_dea_1);
		tx   = dnbinom_mu(1, size, st_v, 0)*k_a;
		k_dec  = plogis(T-1,location,scale,1,0)*b;
		ydot[T-1 + g*T] = tx - k_dec*y[T-1 + g*T] + k_dea*y[T-2 + g*T];
	}
}