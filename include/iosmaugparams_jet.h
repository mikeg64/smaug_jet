

//boundary conditions
//period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
// U upper and L lower boundary type for each direction
#define BCU0 4
#define BCL0 4
#define BCU1 4
#define BCL1 4
#define BCU2 4
#define BCL2 4


//Define parameter block
#define PGAMMA 1.66666667
#define PMU 1.0
#define PETA 0.0
#define PGRAV0 -274.0
#define PGRAV1 0.0
#define PGRAV2 0.0


#define PCMAX 0.02
#define PCOURANT 0.15
#define PRKON 0.0
#define PSODIFON 0.0
#define PMODDTON 0.0
#define PDIVBON 0.0
#define PDIVBFIX 0.0
#define PHYPERDIFMOM 1.0
#define PREADINI 1.0
#define PCFGSAVEFREQUENCY 1

#define PMAXVISCOEF 0.0
#define PCHYP3 0.0
#define PCHYP 0.02
#define PCHYPRHO 0.02
#define PCHYPENERGY 0.02
#define PCHYPB0 0.02
#define PCHYPB1 0.02
#define PCHYPB2 0.02
#define PCHYPMOM0 0.2
#define PCHYPMOM1 0.2
#define PCHYPMOM2 0.2


#define METADDIR "out"
#define METADAUTHOR "MikeG"
#define METADSDATE "Nov 2009"
#define METADPLATFORM "swat"
#define METADDESC "A simple test of SAAS"
#define METADNAME "test1"
#define METADINIFILE "test1.ini"
#define METADLOGFILE "test1.log"
#define METADOUTFILE "test1.out"

#define NI 508
#define NJ 124
#define NK 124

#define DT 0.0001
#define NT 101

#define XMAX 24000000.0
#define YMAX 4000000.0
#define ZMAX 4000000.0

#define XMIN 133333.33
#define YMIN 1953.1
#define ZMIN 1953.1




#define CFGFILE "configs/3D_tubeact_128_128_128_asc_50.ini"
#define CFGOUT "/nobackup/users/cs1mkg/proj/smaug_jet/out/tube2"

real cmax[NDIM];
real courantmax;

int ngi=2;
int ngj=2;
int ngk=2;



//Domain definition
// Define the x domain


//#ifdef USE_SAC
int ni= NI;


real xmax= XMAX;
real xmin= XMIN;

//#endif



// Define the y domain



int nj= NJ ;  //BW test

real ymax= YMAX;
real ymin= YMIN;



#ifdef USE_SAC_3D


int nk = NK;
real zmax= ZMAX;
real zmin= ZMIN;

#endif

real *x, *y, *z;





int step=0;
//real tmax = 200;
real tmax = 0.2;
int steeringenabled=1;
int finishsteering=0;
char configfile[300];

char cfgfile[300]= CFGFILE;

char cfgout[300]= CFGOUT;

Params *d_p, *p;
Meta *metad;


State *d_state, *state;

real dt=DT;


int nt= NT;

real *t;


