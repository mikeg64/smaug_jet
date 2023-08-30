

__device__ __host__
int addsourceterms2_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

   int direction;
  int status=0;

   real xc1,xc2,xc3,r1,r2;
   real xxmax,yymax;
   real dx,dy,dz;
   real aa;
   real s_period;
   real tdep;
   real s_rad1,s_rad2;
   real delta_x, delta_z, delta_t;
   real exp_x,exp_z,exp_xyz;
   real vvv,vvz,AA;

   real n1,n2;
	real qt;





   real xp,yp,zp;
   int i,j,k;
 	  
	  i=ii[0];
	  j=ii[1];
	  k=ii[2];

//printf("source %d %d %d\n",i,j,k);



real xc1Mm=0.5;  // !Mm        z axis
//real xc3Mm=0.5;  // !Mm        z axis
//xc1Mm=1.1;  // !Mm        z axis
//xc1Mm=1.7;  // !Mm        z axis
//real xc2Mm=1.998046;  // !0.99d0  !Mm        x axis
real xc2Mm=1.998046;  // !0.99d0  !Mm        x axis
real xc3Mm=1.998046;  // !Mm        z axis

xc1=xc1Mm*1.0e6;//  !m        z axis
xc2=xc2Mm*1.0e6;//  !m        x axis
xc3=xc3Mm*1.0e6;//  !m        x axis

//yymax=4.0e6;
//xxmax=4.0e6;

//yymax=4.01572e6;
//xxmax=4.01572e6;

yymax=p->xmax[2];
xxmax=p->xmax[1];


xp=(p->xmin[0])+(((real)i)*(p->dx[0]));
 
  
r2=(yp-xc2)*(yp-xc2);

//xxmax=4.0e6;
xxmax=4.01572e6;

exp_x=exp(-r2/(delta_x*delta_x));
delta_x=0.004e6;


qt=p->qt;
r1=(xp-xc1)*(xp-xc1);
delta_z=0.004e6;
n1=2;
n2=2;
//AA=138.89;//*(n1+1)*(n2+1);
//AA=10;//*(n1+1)*(n2+1);
AA=500;//same energy delivered by all drivers
//AA=0;
//AA=0.0356825*(n1+1)*(n2+1);
//AA=2000;// ! 2 km/sec
s_period=300.e0;
yp=(p->xmin[1])+(((real)j)*(p->dx[1]));
zp=(p->xmin[2])+(((real)k)*(p->dx[2]));

exp_z=exp(-r1/(delta_z*delta_z));
exp_xyz=sin(PI*yp*(n1+1)/xxmax)*sin(PI*zp*(n2+1)/yymax)*exp_z;

tdep=sin(qt*2.0*PI/s_period);
vvz=AA*exp_xyz*tdep;


w[fencode3_MODID(p,ii,mom1)]=w[fencode3_MODID(p,ii,mom1)]+(p->dt)*vvz*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
w[fencode3_MODID(p,ii,energy)]=w[fencode3_MODID(p,ii,energy)]+(p->dt)*vvz*vvz*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)])/2.0;
//w[fencode3_MODID(p,ii,mom1)]=w[fencode3_MODID(p,ii,mom1)]+(p->dt)*vvz*w[fencode3_MODID(p,ii,rhob)];

//if(j==63 && i==32  && k==63 /*((p->it)%100 == 0)*/ )
//   printf("Driver amplitude %d %g %g %g %g %g\n",p->it,qt,vvz,w[fencode3_MODID(p,ii,mom1)]/(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]),p->dt,w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);


//if(j==63 && i>=8 && i<=32  && k==63 /*((p->it)%100 == 0)*/ )
   //printf("Driver amplitude %d %g %g %g %g %g\n",p->it,qt,vvz,w[fencode3_MODID(p,ii,mom1)]/(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]),p->dt,w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
//	printf("Driver amplitude %d %g %g %g %g %g\n",p->it,qt,vvz,exp_xyz,p->dt,tdep);




return ( status);
  
  
  
  
}

__device__ __host__
int addsourceterms1_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {


}

