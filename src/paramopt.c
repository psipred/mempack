/***************************************************
 *                    PARAMOPT                     *
 *         By David T. Jones    October 2005       *
 *             UCL Bioinformatics Unit             *
 ***************************************************/

// THIS SOFTWARE MAY ONLY BE USED FOR NON-COMMERCIAL PURPOSES. PLEASE CONTACT
// THE AUTHOR IF YOU REQUIRE A LICENSE FOR COMMERCIAL USE.

#define Version		"0.2"
#define Last_Edit	"16th October 2005"

#define ELITIST

/*
 * This program searches for an optimal set of command line parameters using a
 * genetic-style search.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include "globals.h"
#include "draw_graphs.h"

#define MAXPARAMS 1000

/* Utility definitions */

#define FALSE 0
#define TRUE 1
#define BIG (1000000)
#define VBIG (1e32F)

#define SQR(x) ((x)*(x))
#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))
#define CH alloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define vecadd(a,b,c) (a[0]=b[0]+c[0],a[1]=b[1]+c[1],a[2]=b[2]+c[2])
#define vecsub(a,b,c) (a[0]=b[0]-c[0],a[1]=b[1]-c[1],a[2]=b[2]-c[2])
#define vecscale(a,b,c) ((a[0]=b[0]*(c)),(a[1]=b[1]*(c)),(a[2]=b[2]*(c)))
#define vecprod(a,b,c) ((a[0]=b[1]*c[2]-b[2]*c[1]),(a[1]=b[2]*c[0]-b[0]*c[2]),(a[2]=b[0]*c[1]-b[1]*c[0]))
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])
#define veccopy(a,b) ((a[0]=b[0]),(a[1]=b[1]),(a[2]=b[2]))

double  minparam[MAXPARAMS], maxparam[MAXPARAMS];
double vbest = VBIG;
int ncalls = 0;
short   paramtype[MAXPARAMS];

char    progname[512];


typedef struct
{
    double         *genome;
    float           perfval, selval;
    short           evalflg;
}
Schema;

Schema  *curpool, *newpool;
int      poolsize, genlen, *samparr, besti, rotation;
float    mutrate, crosrate, mutscfac;
double   worst, best, avc_perf;


/* Dump a rude message to standard error and exit */
void fail(char *fmt, ...)
{
    va_list ap;
    
    va_start(ap, fmt) ;
    fprintf(stderr, "*** ");
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    
    exit(-1);
}

/* Marsaglia Universal Double Precision Float RNG */

double U[98];

double uni64()
{
    double x;
    static double c = 0.0;
    static int i=97, j=33;
    
    x = U[i] - U[j];

    if (x < 0.0)
	x += 1.0;

    U[i] = x;
    
    if (! --i)
	i = 97;

    if (! --j)
	j = 97;
    
    c -= 362436069876.0/9007199254740992.0;

    if (c < 0.0)
	c += 9007199254740881.0/9007199254740992.0;
    
    x -= c;

    if (x >= 0.0)
	return x;

    return x + 1.0;
}


void randomise()
{
    double s,t;
    int x = 0,y = 0,i,j;
    struct timeval tv;

    /* Attempt to generate a random state unique to this process/host/time */

    if (!gettimeofday(&tv, NULL))
    {
	x = tv.tv_sec;
	y = tv.tv_usec;
    }
    else
	fail("randomise: cannot generate random number seeds!");

    x ^= y;
    
    x ^= (unsigned int)gethostid();
    y ^= (unsigned int)getpid();

    for (i=1; i<98; i++) {
	s=0.0;
	t=0.5;
	
	for(j=1; j<=53; j++)
	{
	    x=(6969*x)%65543;
	    y=(8888*x)%65579;
	    if (((x^y)&32)>0)
		s += t;
	    t *= 0.5;
	}
	
	U[i]=s;
    }
}


/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((int) ((low) + ((high)-(low)+1) * uni64()))

/* Generate gaussian deviate with mean 0 and stdev 1 */
double gaussrnd()
{
    double x1, x2, w;
 
    do {
	x1 = 2.0 * uni64() - 1.0;
	x2 = 2.0 * uni64() - 1.0;
	w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);

    return x1 * w;
}


/* Evaluate given set of parameter values */
double  eval(double *t)
{
    int i;
    
    all_rotations[rotation].clear();	
    for (i = 0; i < genlen; i++){
    	all_rotations[rotation].push_back((int)t[i]);
    }    

    double v = optimise_rotation(rotation);
	
    ncalls++;
    
    if (v < vbest)
    {
	printf("Best score %f after %d function evaluations.\n", v, ncalls);
	vbest = v;
    }

    return v;
}

/* Initialize the 'world' */
void ga_init(void)
{
    int i, j;
    bool allocd = false;

    /* population arrays */
    if (!allocd)
    {
	curpool = (Schema*) calloc(poolsize, sizeof(Schema));
	newpool = (Schema*) calloc(poolsize, sizeof(Schema));
	samparr = (int*) calloc(poolsize, sizeof(int));

	if (!curpool || !newpool || !samparr)
	    fail("ga_init: cannot create population arrays!");
    }

    for (i = 0; i < poolsize; i++)
    {
	if (!allocd)
	{

	    curpool[i].genome = (double*) calloc(genlen, sizeof(double));
	    newpool[i].genome = (double*) calloc(genlen, sizeof(double));

	    if (curpool[i].genome == NULL || newpool[i].genome == NULL)
		fail("ga_init: cannot create schema!");
	}

	/* Initialize population with random values */
	for (j=0; j<genlen; j++)
	    if (!paramtype[j])
		curpool[i].genome[j] = randint((int)minparam[j], (int)maxparam[j]);
	    else
		curpool[i].genome[j] = uni64() * (maxparam[j] - minparam[j]) + minparam[j];
    
	curpool[i].evalflg = TRUE;
    }
    
    allocd = true;
}

int schcmp(const void *sch1, const void *sch2)
{
    if (((Schema *) sch1)->perfval < ((Schema *) sch2)->perfval)
	return -1;
    else if (((Schema *) sch1)->perfval > ((Schema *) sch2)->perfval)
	return 1;

    return 0;
}

/* Sort pool into ascending order of perfval */
void            sortpool(Schema * pool)
{
    qsort((void *) pool, poolsize, sizeof(Schema), schcmp);
}

/* Select new population from old */
void gaselect(void)
{
    double          ptr;	/* determines fractional selection       */
    double          sum;	/* control for selection loop           */
    double          fitsum;	/* sum of fitness values */
    int             i,k;

#if 1
    sortpool(curpool);

    /* denominator for ordinal selection probabilities */
    for (fitsum = i = 0; i < poolsize; i++)
    {
	curpool[i].selval = poolsize - i;
	fitsum += curpool[i].selval;
	if (curpool[i].perfval == best)
	    besti = i;
    }
#else
    /* denominator for selection probabilities */
    for (fitsum = i = 0; i < poolsize; i++)
    {
	curpool[i].selval = worst - curpool[i].perfval;
	fitsum += curpool[i].selval;
	if (curpool[i].perfval == best)
	    besti = i;
    }
#endif

    for (i = 0; i < poolsize; i++)
    {
	sum = 0.0;
	k = -1;			/* index of next Selected structure */

	ptr = fitsum * uni64();	/* spin the wheel one time */

	do
	{
	    k++;
	    sum += curpool[k].selval;
	}
	while (sum < ptr && k < poolsize - 1);

	samparr[i] = k;
    }

#if 0
    /* randomly shuffle indices to new structures */
    for (i = 0; i < poolsize; i++)
    {
	j = randint(i, poolsize - 1);
	temp = samparr[j];
	samparr[j] = samparr[i];
	samparr[i] = temp;
    }
#endif

    /* Form the new population */
    for (i = 0; i < poolsize; i++)
    {
	k = samparr[i];
	memcpy(newpool[i].genome, curpool[k].genome, genlen * sizeof(double));
	newpool[i].perfval = curpool[k].perfval;
	newpool[i].evalflg = FALSE;
    }

#ifdef ELITIST
    /* Elitist strategy... */
/*    printf("besti = %d %f %f\n", besti, best, curpool[besti].perfval); */
    memcpy(newpool[besti].genome, curpool[besti].genome, genlen * sizeof(double));
    newpool[besti].perfval = curpool[besti].perfval;
    newpool[besti].evalflg = FALSE;
#endif
}


void mutate(float prob)
{
    int i, j;
    double delta;
    static double scale=0.25;

    if (prob > 0.0)
	for (i = 0; i<poolsize; i++)
	    if (i != besti)
		for (j=0; j<genlen; j++)
		    if (uni64() < prob)
		    {
			delta = scale * gaussrnd() * (maxparam[j] - minparam[j]);
			
			if (!paramtype[j] && delta > -1.0 && delta < 1.0 )
			    delta = (delta < 0.0) ? -1.0 : 1.0;

			newpool[i].genome[j] += delta;

			if (newpool[i].genome[j] > maxparam[j])
			    newpool[i].genome[j] = maxparam[j];
			if (newpool[i].genome[j] < minparam[j])
			    newpool[i].genome[j] = minparam[j];

			newpool[i].evalflg = TRUE;
		    }
    /* Apply mutation scaling factor */
    scale *= mutscfac;
}

/* Randomly crossover pool */
void crossovr(void)
{
    int             i,p1;
    double          old1[MAXPARAMS], old2[MAXPARAMS];

    for (i = 0; i < poolsize - 1; i += 2)
	if (uni64() < crosrate)
	{

#ifdef ELITIST
	    if (i == besti || i + 1 == besti)
		continue;
#endif

	    /* Multi point crossover */

	    memcpy(old1, newpool[i].genome, genlen * sizeof(double));
	    memcpy(old2, newpool[i + 1].genome, genlen * sizeof(double));

	    for (p1=0; p1 < genlen; p1++)
		if (uni64() < 0.5)
		{
		    newpool[i].genome[p1] = newpool[i + 1].genome[p1];
		    newpool[i + 1].genome[p1] = old1[p1];
		}

	    if (memcmp(newpool[i].genome, old1, genlen * sizeof(double)))
		newpool[i].evalflg = TRUE;
	    if (memcmp(newpool[i + 1].genome, old2, genlen * sizeof(double)))
		newpool[i + 1].evalflg = TRUE;
	}
}

void statistics(Schema * pool)
{
    int             i;

    for (i = 0; i < poolsize; i++)
	if (pool[i].evalflg)
	{
	    pool[i].perfval = eval(pool[i].genome);
	    pool[i].evalflg = FALSE;
	}

    avc_perf = best = worst = pool[0].perfval;

    for (i = 1; i < poolsize; i++)
    {
	avc_perf += pool[i].perfval;
	if (pool[i].perfval > worst)
	    worst = pool[i].perfval;
	if (pool[i].perfval < best)
	    best = pool[i].perfval;
    }

    avc_perf /= (float) poolsize;
}

void run_ga(void)
{
    int             gen, prevgen = 0,opt_flag;
    float           prevbest;
    Schema         *temp;

    prevbest = VBIG;

    statistics(curpool);
    
    //printf("Initial: %g %g %g\n\n", worst, avc_perf, best);

    for (gen = 1;; gen++)
    {
	gaselect();

	crossovr();
	mutate(mutrate);
	    
	statistics(newpool);
	    
	//printf("%d %g %g %g\n", gen, worst, avc_perf, best);
	fflush(stdout);
	    
	if (best < prevbest)
	{
	    prevbest = best;
	    prevgen = gen;
	    opt_flag = TRUE;
	}

	if (gen - prevgen > 50 || worst == best)
	{
	    puts("*** Convergence detected!");
	    break;
	}
	
	temp = newpool;
	newpool = curpool;
	curpool = temp;
    }
}

/* Read parameters */
void readparams(int helices){

	genlen = helices;
	poolsize = 50;
	mutrate = 0.1;
	crosrate = 0.8;
	mutscfac = 1;

	for (int i=0; i<genlen; i++){
		paramtype[i] = 1;
		minparam[i] = 0;
		maxparam[i] = 359;
   	 }
}

double optimise_parameters(int helices, int r){

    	rotation = r;
    	vbest = VBIG;
	ncalls = 0;

   	 /*
    	printf("Optimum Parameter Search Program\n");
    	printf("Build date : %s\n",__DATE__);
    	printf("Copyright (C) 1994/2005 David T. Jones\n\n");
    	*/

    	randomise();
    	readparams(helices);

    	/*
    	printf("Number of parameters = %d\n", genlen);
    	printf("Pool size = %d\n", poolsize);
    	printf("Mutation rate = %f\n", mutrate);
   	 printf("Crossover rate = %f\n", crosrate);
    	printf("Mutation scaling factor = %f\n", mutscfac);
    	*/

    	ga_init();
    	run_ga();
    	printf("Best score:\t%f\n",best);
	return(best);

}
