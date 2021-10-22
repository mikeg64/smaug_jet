
/*#ifndef HYPERDIF_H_
#define HYPERDIF_H_*/
//	#include <iome/simulation/stdsoap2.h>
    	#include <unistd.h>
	#include <sys/stat.h>
	#include <sys/types.h>
	#include <sys/wait.h>
        #include <sys/time.h>







#include <time.h>
#include "iotypes.h"
//#include "initialisation.h"


#include "initialisation_user.h"






char **hlines;

int ngi,ngj,ngk;
int ni,nj,nk;

real xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz;

real cmax[NDIM];
real courantmax;                
char configfile[300];
char cfggathout[300];
int nt;

struct params *d_p,*p;
struct state *d_state, *state;
struct meta *meta;



// Define time-domain
real dt;


real **d_gw, **d_gwnew;
real **d_gwtemp,**d_gwtemp1,**d_gwtemp2;
real **d_gwmod,  **d_gdwn1,  **d_gwd;


real *d_w;
real *d_wnew;

real *d_wmod,  *d_dwn1,  *d_dwn2,  *d_dwn3,  *d_dwn4,  *d_wd;

real *w,*wnew,*wd,*wdnew, *temp2,*wmod;
real *d_wtemp,*d_wtemp1,*d_wtemp2;

#ifdef USE_MULTIGPU
  real *gmpivisc0,*gmpivisc1,*gmpivisc2, *gmpiw, *gmpiwmod, *gmpiw0, *gmpiwmod0, *gmpiw1, *gmpiwmod1, *gmpiw2, *gmpiwmod2;
  real *d_gmpivisc0,*d_gmpivisc1,*d_gmpivisc2, *d_gmpiw, *d_gmpiwmod, *d_gmpiw0, *d_gmpiwmod0, *d_gmpiw1, *d_gmpiwmod1, *d_gmpiw2, *d_gmpiwmod2;

  real **d_ggmpivisc0,**d_ggmpivisc1,**d_ggmpivisc2, **d_ggmpiw, **d_ggmpiwmod, **d_ggmpiw0, **d_ggmpiwmod0, **d_ggmpiw1, **d_ggmpiwmod1, **d_ggmpiw2, **d_ggmpiwmod2;

  real **d_gmpiviscr0,**d_gmpiviscr1,**d_gmpiviscr2, **d_gmpiwr0, **d_gmpiwmodr0, **d_gmpiwr1, **d_gmpiwmodr1, **d_gmpiwr2, **d_gmpiwmodr2, **d_gmpiwmodr2;


#endif


#ifdef USE_MULTIGPU
//buffers to use on GPU
  real *d_gmpisendbuffer;
  real *d_gmpirecvbuffer;

   
  real *d_gmpisrcbufferl;
  real *d_gmpisrcbufferr;
  real *d_gmpitgtbufferl;
  real *d_gmpitgtbufferr;
#endif





/*----------------------*/ 
real second()
{

   /*REAL secs;
   clock_t Time;
   Time = clock();

   secs = (real) Time / (real) CLOCKS_PER_SEC;
   return secs;*/
   real retval;
	static long zsec=0;
	static long zusec=0;
	real esec;
	
	struct timeval tp;
	//struct timezone tzp;
	
	gettimeofday(&tp, NULL);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	return retval;

}

int fencode3_test (struct params *dp,int *ii, int field);
int encode3_in(struct params *dp,int ix, int iy, int iz, int field);
void initconfig(struct params *k, struct meta *md, real *w, real *wd);

#ifdef USE_MPI
	int sacencodempivisc0 (struct params *p,int ix, int iy, int iz, int bound,int dim);
	int sacencodempivisc1 (struct params *p,int ix, int iy, int iz, int bound,int dim);
	int sacencodempivisc2 (struct params *p,int ix, int iy, int iz, int bound,int dim);
	int sacencodempiw0 (struct params *p,int ix, int iy, int iz, int field,int bound);
	int sacencodempiw1 (struct params *p,int ix, int iy, int iz, int field,int bound);
	int sacencodempiw2 (struct params *p,int ix, int iy, int iz, int field,int bound);
	int encode3p2_sacmpi (struct params *dp,int ix, int iy, int iz, int field);

	void mgpuinit(struct params *p);
	void mgpufinalize(struct params *p);

	void mgpusetnpediped(struct params *p, char *string);
	void ipe2iped(struct params *p);

	void iped2ipe(int *tpipe,int *tpnp, int *oipe);

	void mgpuneighbours(int dir, struct params *p);

	void mpisend(int nvar,real *var, int *ixmin, int *ixmax  ,int qipe,int iside, int dim, struct params *p);

	void mpisendmod(int nvar,real *var, int *ixmin, int *ixmax  ,int qipe,int iside, int dim, struct params *p);
	void mpirecvbuffer(int nvar,int *ixmin, int *ixmax  ,int qipe,int iside,int dim, struct params *p);
	void mpirecvbuffermod(int nvar,int *ixmin, int *ixmax  ,int qipe,int iside,int dim, struct params *p);
	void mpibuffer2var(int iside,int nvar,real *var, int *ixmin, int *ixmax, int dim, struct params *p);

	void mpibuffer2varmod(int iside,int nvar,real *var, int *ixmin, int *ixmax, int dim, struct params *p);

	void gpusync();

	void mpibound(int nvar,  real *var1, real *var2, real *var3, struct params *p, int idir);

	void mpiboundmod(int nvar,  real *var1, real *var2, real *var3, struct params *p, int idir);

	void mpireduce(real *a, void * mpifunc);  //nb void* should be MPI_Op

	void mpiallreduce(real *a, void * mpifunc);

	void mpivisc( int idim,struct params *p, real *var1, real *var2, real *var3);

#endif








