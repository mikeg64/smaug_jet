













ngi=2;
ngj=2;
ngk=2;


//Domain definition
// Define the x domain
#ifdef USE_SAC
//vac ozt
//int ni;
ni=252+2*ngi;    //OZT tests
//ni+=2*ngi;
//ni=512;
//real xmax = 6.2831853;  
xmax=1.0;
dx = xmax/(ni);
#endif



// Define the y domain
#ifdef USE_SAC
//vac ozt
nj = 252+2*ngj;  //OZT tests
//int nj=2;  //BW test
//nj+=2*ngj;
//nj=512;
//real ymax = 6.2831853; 
ymax = 1.0;   
dy = ymax/(nj);    
//nj=41;
#endif






  



int step=0;
//real tmax = 200;
real tmax = 0.2;
int steeringenabled=1;
int finishsteering=0;

//char *cfgfile="zero1.ini";

//char *cfgfile="zero1_np020203.ini";
//char *cfgfile="zero1_np0201.ini";
//char *cfgfile="configs/zero1_ot_asc.ini";
char *cfgfile="zero1_ot_asc_256.ini";
//char *cfgfile="zero1_BW_bin.ini";
//char *cfgout="zero1_np010203."
char *cfgout="out/zeroOT";
//char *cfgout="zero1_np0201.out";






#ifdef USE_SAC
dt=0.0002;  //OZT test
#endif


//nt=3000;
//nt=5000;
//nt=200000;
//nt=150000;
nt=101;


real *t=(real *)calloc(nt,sizeof(real));
printf("runsim 1%d \n",nt);
//t = [0:dt:tdomain];
for(i=0;i<nt;i++)
		t[i]=i*dt;

//real courant = wavespeed*dt/dx;

p->n[0]=ni;
p->n[1]=nj;
p->ng[0]=ngi;
p->ng[1]=ngj;




p->dt=dt;
p->dx[0]=dx;
p->dx[1]=dy;







//ozt test
p->gamma=1.66667;  //OZ test
//p->gamma=2.0;  //BW test
//p->gamma=5.0/3.0;  //BACH3D
//alfven test
//p->gamma=1.4;






p->mu=1.0;
p->eta=0.0;
p->g[0]=0.0;
p->g[1]=0.0;
p->g[2]=0.0;
#ifdef USE_SAC_3D

#endif
//p->cmax=1.0;
p->cmax=0.02;

p->rkon=0.0;
p->sodifon=0.0;
p->moddton=0.0;
p->divbon=0.0;
p->divbfix=0.0;
p->hyperdifmom=1.0;
p->readini=1.0;
p->cfgsavefrequency=1;


p->xmax[0]=xmax;
p->xmax[1]=ymax;
p->nt=nt;
p->tmax=tmax;
p->steeringenabled=steeringenabled;
p->finishsteering=finishsteering;

p->maxviscoef=0;
//p->chyp=0.0;       
//p->chyp=0.00000;
p->chyp3=0.00000;


for(i=0;i<NVAR;i++)
  p->chyp[i]=0.0;

p->chyp[rho]=0.2;
p->chyp[energy]=0.2;
p->chyp[b1]=0.2;
p->chyp[b2]=0.2;
p->chyp[mom1]=0.2;
p->chyp[mom2]=0.2;
p->chyp[rho]=0.2;

p->npe=1;


#ifdef USE_MPI
//number of procs in each dim mpi only
p->pnpe[0]=2;
p->pnpe[1]=1;
p->pnpe[2]=1;
#endif


//iome elist;
//meta meta;







//elist.server=(char *)calloc(500,sizeof(char));




	//strcpy(elist.server,"localhost1");
	//elist.port=80801;
	//elist.id=0;


