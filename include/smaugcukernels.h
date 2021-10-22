
#include "iotypes.h"
//#include "iobparams.h"




extern "C" int cuinit(struct params **p, struct bparams **bp, real **wmod,real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
extern "C" int cuupdatemod(struct params **p, struct bparams **bp,real **w, real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
extern "C" int initgrid(struct params **p,   struct state **state, real **wd, struct params **d_p, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

extern "C" int cufinish(struct params **p, real **w, real **wnew,   struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

extern "C" int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp,struct state **d_s,  real **d_wmod,  int order,int idir,int field);
extern "C" int cuupdate(struct params **p, real **w, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_w,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);

extern "C" int cuupdatehostwd(struct params **p, real **wd, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_wd,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);
extern "C" int cuupdatedevicewd(struct params **p, real **wd, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_wd,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);



extern "C" int cucentdiff1(struct params **p, struct params **d_p, struct state **d_s,real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field, int dir);
extern "C" int cucentdiff2(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt, int field,int dir);

extern "C" int cugrav(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);
extern "C" int cusource(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);

extern "C" int cuadvance(struct params **p, struct params **d_p,    real **d_wmod, real **d_w,int order);
extern "C" int cucomputedervfields(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order);
extern "C" int cucomputevels(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
extern "C" int cucomputemaxc(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir, real **wd, real **d_wtemp);
extern "C" int cucomputemaxcourant(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir, real **wd, real **d_wtemp);
extern "C" int cucomputec(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
extern "C" int cucomputept(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
extern "C" int cucomputepk(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
extern "C" int cucomputepbg(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
extern "C" int cudivb(struct params **p,struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);

extern "C" int cuhyperdifmomsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt);

extern "C" int cuhyperdifmomsourcene1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt);

extern "C" int cuhyperdifesource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real **d_wtemp, int field, int dim, real dt);

extern "C" int cuhyperdifbsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt);

extern "C" int cuhyperdifbsourcene1(struct params **p, struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt);


extern "C" int cuhyperdifrhosource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, real dt);



extern "C" int cunushk1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
extern "C" int cugetdtvisc1(struct params **p,  struct params **d_p,   real **d_wmod, real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

extern "C" int cuhyperdifvisc1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand);
extern "C" int cuhyperdifvisc1r(struct params **p,  struct params **d_p,   real **d_wmod,real **wd,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
extern "C" int cuhyperdifvisc1l(struct params **p,  struct params **d_p,   real **d_wmod, real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
extern "C" int cuhyperdifvisc1ir(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
extern "C" int cuhyperdifvisc1il(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);


extern "C" int cuhyperdifviscmax(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, int field, int dim);
extern "C" int cusync(struct params **p);
#ifdef USE_MULTIGPU
 extern "C" int cuinitmgpubuffers(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2);
 extern "C" int cufinishmgpu(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2);
extern "C" int cucopywtompiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);
extern "C" int cucopywfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);

extern "C" int cucopywdtompiwd(struct params **p,real **wd,     real **gmpiw0,    real **gmpiw1,    real **gmpiw2, struct params **d_p  ,real **d_wd,    real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2,  int order,int idir);

extern "C" int cucopywdfrommpiwd(struct params **p,real **wd,     real **gmpiw0,    real **gmpiw1,    real **gmpiw2, struct params **d_p  ,real **d_wd,    real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2,  int order, int idir);


extern "C" int cucopytompivisc(struct params **p,real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2);
extern "C" int cucopyfrommpivisc(struct params **p,real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2);
extern "C" int cusetgpu(struct params **p);



extern "C" int cucopywtompiwmod(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);

extern "C" int cucopywmodfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);


#endif
