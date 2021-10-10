

__device__ __host__
int addsourceterms2_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;

   real xc1,xc2,r1,r2;
   real xxmax,yymax;
   real dx,dy,dz;
   real aa;
   real s_period;
   real qt, tdep;
   real s_rad1,s_rad2;
   real exp_z, exp_x, exp_xy, exp_y;

   real vvv;

   real tim, timtanh;



   real xp,yp,zp;
   int i,j,k;
   int n1;
 	  
	  i=ii[0];
	  j=ii[1];

     xc1=4.0e6;
    xc2=300000.0;
    qt=p->qt;

    //if(qt<300)
    //   xc2=xc2+125000;

    aa=1000;
    n1=1;

          xp=(p->xmin[0])+(((real)i)*(p->dx[0]));
          yp=(p->xmin[1])+(((real)j)*(p->dx[1]));
     
   xxmax=(p->xmax[0])-2*((p->dx[0]));
    s_period=100.0;
    tdep=exp(-(qt-s_period)*(qt-s_period));
    //tdep=sin(qt*2.0*PI/s_period);

     r1=(xp-xc1)*(xp-xc1);
     r2=(yp-xc2)*(yp-xc2);

 
 s_rad1=1870.0; // from Mackenzie dover app 913:19(10pp), 2021 May 20
 s_rad2=1000000.0;
 exp_x=exp(-(r1/(s_rad1*s_rad1)));
 exp_y=exp(-(r2/(s_rad2*s_rad2)));
 //exp_xy=sin(PI*xp*(n1+1)/xxmax)*exp_z;

   tim=PI+(PI*((qt-s_period)/s_period)); 
   timtanh=(exp(2*tim)-1)/(exp(2*tim)+1);
   exp_xy=-(1+timtanh)*exp_x*exp_y/2;
   

      vvv=aa*tdep*exp_xy;

     // if(i==3 && j==149)
      //if(i==3 && j==200)
     /* {
                  p->test=vv;
               p->chyp[0]=xp;
                p->chyp[1]=yp;
       }*/

       /* if(i==9 && j==63 && k==63) 
	{
                p->test=(w[fencode3_MODID(p,ii,rho)]);
                p->chyp[0]=vx;
                p->chyp[1]=vy;
                p->chyp[2]=(w[fencode3_MODID(p,ii,mom1)]);
	}*/


//if(i==512 && j==15 )
//printf("%g %g  %g %g %g \n", xp, yp,vvv,r1, r2);
 

// if(j==2 || j==3)
//{
                           w[fencode3_MODID(p,ii,mom2)]=w[fencode3_MODID(p,ii,mom2)]+(p->dt)*vvv*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
  
                          w[fencode3_MODID(p,ii,energy)]=w[fencode3_MODID(p,ii,energy)]+(p->dt)*(vvv*vvv)*(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)])/2.0;
//}
  

  return ( status);
}

__device__ __host__
int addsourceterms1_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

   

}

