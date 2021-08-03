
__device__ __host__
real usersource2a_cd1(real *dw, real *wd, real *w, struct params *p,int *ii, int dir) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real source=0;


        #if defined USE_SAC  || defined USE_SAC_3D
//computept3_cd1(w,wd,p,ii);



  /*    		source= -wd[fencode3_cd1(p,ii,ptb)]*grad3d_cd1(wd,p,ii,vel1+dir,dir);
               source += +w[fencode3_cd1(p,ii,b1b)]*w[fencode3_cd1(p,ii,b1b+dir)]*grad3d_cd1(wd,p,ii,vel1,0)+w[fencode3_cd1(p,ii,b2b)]*w[fencode3_cd1(p,ii,b1b+dir)]*grad3d_cd1(wd,p,ii,vel1+1,1);*/
         #endif


        #if defined USE_SAC_3D
            //   source += +w[fencode3_cd1(p,ii,b3b)]*w[fencode3_cd1(p,ii,b1b+dir)]*grad3d_cd1(wd,p,ii,vel3,2);
        #endif

  return source;


 
}

__device__ __host__
real usersource1a_cd1(real *dw, real *wd, real *w, struct params *p,int *ii, int dir) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real source=0;


        #if defined USE_SAC  || defined USE_SAC_3D
//computept3_cd1(w,wd,p,ii);



      	/*	source= -wd[fencode3_cd1(p,ii,ptb)]*grad3d_cd1(wd,p,ii,vel1+dir,dir);
               source += +w[fencode3_cd1(p,ii,b1b)]*w[fencode3_cd1(p,ii,b1b+dir)]*grad3d_cd1(wd,p,ii,vel1,0)+w[fencode3_cd1(p,ii,b2b)]*w[fencode3_cd1(p,ii,b1b+dir)]*grad3d_cd1(wd,p,ii,vel1+1,1);*/
         #endif


        #if defined USE_SAC_3D
              // source += +w[fencode3_cd1(p,ii,b3b)]*w[fencode3_cd1(p,ii,b1b+dir)]*grad3d_cd1(wd,p,ii,vel3,2);
        #endif

  return source;


}



__device__ __host__
int addsourceterms2_cd1(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divsource=0;
  //dw[fencode3_cd1(p,ii,field)]= grad_cd1(wd,p,ii,source,dir);//+grad_cd1(wd,p,ii,f2,1); 


 #if defined USE_SAC  ||  defined USE_SAC_3D

  
 /* if(field==energy)
  {    
     computept3_cd1(w,wd,p,ii);
     dw[fencode3_cd1(p,ii,field)]=usersource2a_cd1(dw, wd, w, p,ii,dir)+w[fencode3_cd1(p,ii,rho)]*((p->g[dir])*w[fencode3_cd1(p,ii,mom1+dir)]    )/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]);
   }*/


 #endif
  return ( status);
}

__device__ __host__
int addsourceterms1_cd1(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divsource=0;
  //dw[fencode3_cd1(p,ii,field)]= grad_cd1(wd,p,ii,source,dir);//+grad_cd1(wd,p,ii,f2,1); 


 #if defined USE_SAC  ||  defined USE_SAC_3D

  
  /*if(field==energy)
  {    
     computept3_cd1(w,wd,p,ii);
     dw[fencode3_cd1(p,ii,field)]=usersource2a_cd1(dw, wd, w, p,ii,dir)+w[fencode3_cd1(p,ii,rho)]*((p->g[dir])*w[fencode3_cd1(p,ii,mom1+dir)]    )/(w[fencode3_cd1(p,ii,rho)]+w[fencode3_cd1(p,ii,rhob)]);
   }*/


 #endif
  return ( status);
}

