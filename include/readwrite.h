#include "iotypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int createlog(char *logfile);
int appendlog(char *logfile, struct params p, struct state s);
int writeconfig(char *name,int n,struct params p, struct meta md, real *w);
int writevtkconfig(char *name,int n,struct params p, struct meta md, real *w);
int writevacconfig(char *name,int n,struct params p, struct meta md, real *w,real *wd, struct state st);
int writevacgatherconfig(char *name,int n,struct params p, struct meta md, real *w,real *wd, struct state st);
int readconfig(char *cfgfile, struct params p, struct meta md, real *w);
int readasciivacconfig(char *cfgfile, struct params p, struct meta md, struct state *st, real *w,real *wd, char **hlines, int mode);
/*Big problems with reading fortran unformatted "binary files" need to include 
  record field*/
int readbinvacconfig(char *name,struct params p, struct meta md, real *w,real *wd, struct state st);
int writeasciivacconfig(char *cfgfile, struct params p, struct meta md, real *w,real *wd, char **hlines, struct state st, int mode);
int createconfigsegment(struct params p,  real *wnew,real *wdnew, real *w,real *wd);
int gathersegment(struct params p,  real *wnew,real *wdnew, real *w,real *wd);
void readatmos(struct params p,real *w);
