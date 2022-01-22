// iosmaug.cpp : Main routine for GPU enabled SAC
#include "../include/iosmaug.h"
//#include "../include/readwrite.h"
//#include "../include/initialisation.h"
//#include "../include/iotypes.h"


/*#include "../include/defs.h"
#include "../include/iobparams.h"*/






#ifdef USE_MPI
	#include "mpi.h"

MPI_Comm *comm;
MPI_Request *request;
double gwall_time;
real *gmpisendbuffer;
real *gmpirecvbuffer;


real **gmpisrcbufferl;
real **gmpisrcbufferr;
real **gmpitgtbufferl;
real **gmpitgtbufferr;


int gnmpirequest,gnmpibuffer,gnmpibuffermod;
MPI_Request *gmpirequest;



#endif



int main(int argc, char* argv[])
{



int itype=-1;
int status=1;
int mode=run;//run a model 1=scatter 2=gather
int it=0; //test integer to be returned 
int n;
//getintparam_( int elist.id,char *sname,int *iv,  int elist.port, char *selist.server );
//int elist.id=0;
//int elist.port=8080;

int i1,i2,i3,j1;
int i,j,k,iv;

	char ext[4];
	char tcfg[300];
	char stemp[300];


	char *pch1,*pch2;



char *sdir=(char *)calloc(500,sizeof(char));
char *name=(char *)calloc(500,sizeof(char));
char *outfile=(char *)calloc(500,sizeof(char));
char *formfile=(char *)calloc(500,sizeof(char));
hlines = (char **) malloc(5*sizeof(char*));

char configinfile[300];
char cfggathout[300];

 real tcom,tcom1, tcom2,tv,tcal,tc;

tcom=0.0;
tcal=0.0;




struct bparams *d_bp;
struct bparams *bp=(struct bparams *)malloc(sizeof(struct bparams));

p=(struct params *)malloc(sizeof(struct params));
state=(struct state *)malloc(sizeof(struct state));
meta=(struct meta *)malloc(sizeof(struct meta));


FILE *portf;

#include "../include/iosmaugparams.h"


if(argc>3  && strcmp(argv[2],"gather")==0 && (atoi(argv[3])>=0) && (atoi(argv[3])<=nt)) 
{    
  mode=gather;
}

if(argc>2  && strcmp(argv[2],"scatter")==0)
{
  mode=scatter;
}

if(argc>2  && strcmp(argv[2],"init")==0)
{
  mode=init;
  printf("init mode=3\n");
}
p->mode=mode;




       /*********************************************************************************************************/
       /* Start of section to set domain sizes and config filenames*/
       /*********************************************************************************************************/

real *x=(real *)calloc(ni,sizeof(real));
for(i=0;i<ni;i++)
		x[i]=i*dx;

real *y=(real *)calloc(nj,sizeof(real));
for(i=0;i<nj;i++)
		y[i]=i*dy;


//set boundary types
for(int ii=0; ii<NVAR; ii++)
for(int idir=0; idir<NDIM; idir++)
for(int ibound=0; ibound<2; ibound++)
{
   (p->boundtype[ii][idir][ibound])=0;  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
}


/*meta.directory=(char *)calloc(500,sizeof(char));
meta.author=(char *)calloc(500,sizeof(char));
meta.sdate=(char *)calloc(500,sizeof(char));
meta.platform=(char *)calloc(500,sizeof(char));
meta.desc=(char *)calloc(500,sizeof(char));
meta.name=(char *)calloc(500,sizeof(char));
meta.ini_file=(char *)calloc(500,sizeof(char));
meta.log_file=(char *)calloc(500,sizeof(char));
meta.out_file=(char *)calloc(500,sizeof(char));

strcpy(meta.directory,"out");
strcpy(meta.author,"MikeG");
strcpy(meta.sdate,"Nov 2009");
strcpy(meta.platform,"swat");
strcpy(meta.desc,"A simple test of SAAS");
strcpy(meta.name,"test1");
strcpy(meta.ini_file,"test1.ini");
strcpy(meta.log_file,"test1.log");
strcpy(meta.out_file,"test1.out");*/



	strcpy(stemp,cfgfile);


	pch1 = strtok (stemp,".");


	sprintf(tcfg,"%s",pch1);
	pch2 = strtok (NULL,".");

	sprintf(ext,"%s",pch2);



	sprintf(configfile,"%s",cfgout);


	#ifdef USE_MULTIGPU
	#ifdef USE_MPI
	     MPI_Init(&argc, &argv);
	#endif

	mgpuinit(p);
	//mgpuinit_stage1(p);



	if(mode==run)
	{
		ipe2iped(p); 

    
		mgpuneighbours(0,p);
		mgpuneighbours(1,p);
	}


	//compute the max and min domain dimensions for each processor
	p->xmax[0]=xmin+(1+(p->pipe[0]))*(xmax-xmin)/(p->pnpe[0]);
	p->xmax[1]=ymin+(1+(p->pipe[1]))*(ymax-ymin)/(p->pnpe[1]);
	p->xmin[0]=xmin+(p->pipe[0])*(xmax-xmin)/(p->pnpe[0]);
	p->xmin[1]=ymin+(p->pipe[1])*(ymax-ymin)/(p->pnpe[1]);

        p->dx[0]=(p->xmax[0]-p->xmin[0])/(p->n[0]);
	p->dx[1]=(p->xmax[1]-p->xmin[1])/(p->n[1]);


	//store global values for max and min domain dimensions
	p->gxmax[0]=xmax;
	p->gxmin[0]=xmin;
	p->gxmax[1]=ymax;
	p->gxmin[1]=ymin;

	#ifdef USE_SAC_3D
	mgpuneighbours(2,p);
	p->xmax[2]=zmin+(1+(p->pipe[2]))*(zmax-zmin)/(p->pnpe[2]);
	p->xmin[2]=zmin+(p->pipe[2])*(zmax-zmin)/(p->pnpe[2]);
	p->gxmax[2]=zmax;
	p->gxmin[2]=zmin;
        p->dx[2]=(p->xmax[2]-p->xmin[2])/(p->n[2]);
	#endif



	  sprintf(configinfile,"%s",cfgfile);

	//adopt the sac MPI naming convention append the file name npXXYY where XX and YY are the
	//number of processors in the x and y directions
	#ifdef USE_MPI
	     #ifdef USE_SAC_3D

		      if(p->ipe>99)
			sprintf(configinfile,"%s_np0%d0%d0%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
		      else if(p->ipe>9)
			sprintf(configinfile,"%s_np0%d0%d0%d_0%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);
		      else
			sprintf(configinfile,"%s_np0%d0%d0%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe,ext);  	     
	     #else
		      if(p->ipe>99)
			sprintf(configinfile,"%s_np0%d0%d_%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
		      else if(p->ipe>9)
			sprintf(configinfile,"%s_np0%d0%d_0%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);
		      else
			sprintf(configinfile,"%s_np0%d0%d_00%d.%s",tcfg,p->pnpe[0],p->pnpe[1],p->ipe,ext);  	     	     
	     #endif
	#endif

	//if doing a scatter or gather set the domain size correctly
	//take a distribution and distribute domain to processors
	if(mode==scatter )
	{
	  printf("Scatter %s %d %d %d\n",cfgfile,p->pnpe[0],p->pnpe[1],p->pnpe[2]);
	  sprintf(configinfile,"%s",cfgfile);
	  p->n[0]=ni*(p->pnpe[0]);
	  p->n[1]=nj*(p->pnpe[1]);
	   #ifdef USE_SAC_3D
		    p->n[2]=nk*(p->pnpe[2]);
	   #endif
	}

	if( mode==gather)
	{
	   ni=ni*(p->pnpe[0]);
	   nj=nj*(p->pnpe[1]);
	   #ifdef USE_SAC_3D
		   nk=nk*(p->pnpe[2]);
	   #endif
	}


	if(mode==init)
	{
	    p->n[0]=ni;
	    p->n[1]=nj;
	    #ifdef USE_SAC_3D
	      p->n[2]=nk;
	    #endif
	}
	printf("config files\n%s \n %s %d %d\n",configinfile,configfile,p->n[0],p->n[1]);


	#else
	     sprintf(configinfile,"%s",cfgfile);
	#endif   //#ifdef USE_MULTIGPU

char *method=NULL;


       /*********************************************************************************************************/
       /* Start of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/
    #ifdef USE_IOME
        if(argc>2)
        {
          //simfile already read by 
          readsim(p,&meta,argv[2],elist);
          //if((p->readini)!=0)
          //   readconfig(meta.ini_file,*p,meta,w);
        }
        else
	  createsim(*p,meta,simfile,elist);

	sprintf(simfile,"%s.xml",meta.name);
        sprintf(newsimfile,"%s_update.xml",meta.name);
     #endif
       /*********************************************************************************************************/
       /* End of section initialising steering and auto metadata collection*/
       /*********************************************************************************************************/


       /*********************************************************************************************************/
       /* Start of section creating arrays on the host*/
       /*********************************************************************************************************/
  	printf("Creating arrays on the host\n");

       #ifdef USE_MULTIGPU
       if(mode==0)
       {
		if( ((p->pnpe[0])>1) &&  (p->pipe[0])==0) (p->n[0])+=ngi;
		if( ((p->pnpe[0])>1) &&  (p->pipe[0])==((p->pnpe[0])-1)) (p->n[0])+=ngi;
		if( ((p->pnpe[1])>1) &&  (p->pipe[1])==0) (p->n[1])+=ngj;
		if( ((p->pnpe[1])>1) &&  (p->pipe[1])==((p->pnpe[1])-1)) (p->n[1])+=ngj;

		#ifdef USE_SAC_3D
			if( ((p->pnpe[2])>1) &&  (p->pipe[2])==0) (p->n[2])+=ngk;
			if( ((p->pnpe[2])>1) &&  (p->pipe[2])==((p->pnpe[2])-1)) (p->n[2])+=ngk;
		#endif
	}
       #endif

	//allocate arrays to store fields, updated fields, dervived quantities and updated derived qunatities
	#ifdef USE_SAC_3D
		wnew=(real *)calloc(ni*nj*nk*NVAR,sizeof(real ));

		wdnew=(real *)calloc(ni*nj*nk*NDERV,sizeof(real ));
		wd=(real *)calloc(((p)->n[0])*((p)->n[1])*((p)->n[2])*NDERV,sizeof(real ));
		wmod=(real *)calloc(2*(1+(((p)->rkon)==1))*((p)->n[0])*((p)->n[1])*((p)->n[2])*NVAR,sizeof(real ));
	#else
		wnew=(real *)calloc(ni*nj*NVAR,sizeof(real ));

		wdnew=(real *)calloc(ni*nj*NDERV,sizeof(real ));
		wd=(real *)calloc(((p)->n[0])*((p)->n[1])*NDERV,sizeof(real ));
		wmod=(real *)calloc(2*(1+(((p)->rkon)==1))*((p)->n[0])*((p)->n[1])*NVAR,sizeof(real ));
	#endif

        #ifdef USE_MULTIGPU
          //parameters used to set sizes of MPI communications buffers and
	  //data storage areas
          int szw,szw0,szw1,szw2,szvisc0,szvisc1,szvisc2;
	  #ifdef USE_SAC
		  szw=4*(  ((p)->n[1])  +  ((p)->n[0])   );
		  szw0=4*NDERV*(  ((p)->n[1])     );
		  szw1=4*NDERV*(  ((p)->n[0])     );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2 )   );
		  szvisc1=4*NVAR*(    (((p)->n[0]) +2 )  );
	  #endif
	  #ifdef USE_SAC_3D	  
		  szw=4*NDERV*(  ((p)->n[1])*((p)->n[2])  +  ((p)->n[0])*((p)->n[2])  +  ((p)->n[0])*((p)->n[1])  );
		  szw0=4*NDERV*(  ((p)->n[1])*((p)->n[2])    );
		  szw1=4*NDERV*(    ((p)->n[0])*((p)->n[2])   );
		  szw2=4*NDERV*(    ((p)->n[0])*((p)->n[1])  );

		  szvisc0=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)  ); 
		  szvisc1=4*NVAR*(   (((p)->n[0])+2)*(((p)->n[2])+2)    );    
		  szvisc2=4*NVAR*(  (((p)->n[1])+2)*(((p)->n[2])+2)   );   
	  #endif




	  #ifdef USE_SAC
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2),sizeof(real));
	  #endif
	  #ifdef USE_SAC_3D
	  temp2=(real *)calloc(NTEMP2*(((p)->n[0])+2)* (((p)->n[1])+2)* (((p)->n[2])+2),sizeof(real));
	  #endif

          //Data areas to store values communicated using MPI
	  gmpiwmod0=(real *)calloc(szw0,sizeof(real));
	  gmpiw0=(real *)calloc(szw0,sizeof(real));
	  gmpiwmod1=(real *)calloc(szw1,sizeof(real));
	  gmpiw1=(real *)calloc(szw1,sizeof(real));

          gmpivisc0=(real *)calloc(szvisc0,sizeof(real));
          gmpivisc1=(real *)calloc(szvisc1,sizeof(real));

	  #ifdef USE_SAC_3D
		  gmpiwmod2=(real *)calloc(szw2,sizeof(real));
		  gmpiw2=(real *)calloc(szw2,sizeof(real));
                  gmpivisc2=(real *)calloc(szvisc2,sizeof(real));
	  #endif
        #endif
       /*********************************************************************************************************/
       /* End of section creating arrays on the host*/
       /*********************************************************************************************************/

	//set initial time step to a large value
	if(p->moddton==1.0)
	{
		p->dt=1.0e-8;
	}
       int its=p->it;


       /*********************************************************************************************************/
       /* Start of section initialising the configuration 
          on the host and on GPU host memory*/
       /*********************************************************************************************************/

       if(mode !=init)
       {
               if((p->readini)==0)
               {
                 printf("init config\n");
		 initconfig(p, &meta, wmod,wd);
                }
		else
                {
	         printf("reading configuration from %s\n",configinfile);
		 readasciivacconfig(configinfile,*p,meta, state,wmod,wd,hlines,mode);
                }
       }


        /*********************************************************************************************************/
       /* Start of section to scatter data
        /*********************************************************************************************************/
        #ifdef USE_MULTIGPU
        //scatter/distribute configuration across each CPU
        if(mode==scatter)
        {
	       gpusync();
               if(p->ipe==0) //currently only processor zero
	       {
                  p->npe=(p->pnpe[0])*(p->pnpe[1])*(p->pnpe[2]);
		  for(i=0; i<p->npe; i++)
		  {
		    p->ipe=i;
                    ipe2iped(p);
		    
		    //copy segment
		    printf("copy segment %d %d %d\n",i,p->npe,p->ipe);                    
		    createconfigsegment(*p, wnew,wdnew,wmod,wd);  //in readwrite.c

		    //writeas
                    //set domain size to size for each processor		   
                    p->n[0]=ni;
                    p->n[1]=nj;
                    #ifdef USE_SAC_3D
                      p->n[2]=nk;
                    #endif
		    writeasciivacconfig(configinfile, *p, meta,  wnew,wdnew, hlines, *state,mode);
                    //set domain size to the global domain size
                    //this will be used when we extract a segment
                    p->n[0]=ni*(p->pnpe[0]);
                    p->n[1]=nj*(p->pnpe[1]);
                    #ifdef USE_SAC_3D
                      p->n[2]=nk*(p->pnpe[2]);
                    #endif
		  }
		}
                gpusync();
        }
       /*********************************************************************************************************/
       /* End of section to scatter data
        /*********************************************************************************************************/
 
/*********************************************************************************************************/
	/* Start of section to gather data
	/*********************************************************************************************************/
	//gather configuration to single output file

	if(mode==gather)
	{
		n=atoi(argv[3]);
		if(p->ipe==0)
		{
			int myipe=p->ipe;
			for(i=0; i<p->npe; i++)
			{
				printf(" here nt=%d pid=%d i=%d\n",n,p->ipe,i);

				p->ipe=i;
				ipe2iped(p);
				strcpy(stemp,cfgout);
				pch1 = strtok (stemp,".");
				sprintf(tcfg,"%s",pch1);

				#ifdef USE_SAC_3D
					if(p->ipe>99)
						sprintf(configinfile,"%s%d_np0%d0%d0%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);
					else if(p->ipe>9)
						sprintf(configinfile,"%s%d_np0%d0%d0%d_0%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);
					else
						sprintf(configinfile,"%s%d_np00%d00%d00%d_00%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->pnpe[2],p->ipe);  	     
				#else
					if(p->ipe>99)
						sprintf(configinfile,"%s%d_np0%d0%d_%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);
					else if(p->ipe>9)
						sprintf(configinfile,"%s%d_np0%d0%d_0%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);
					else
						sprintf(configinfile,"%s%d_np0%d0%d_00%d.out",tcfg,n,p->pnpe[0],p->pnpe[1],p->ipe);  	     	     
				#endif

				//copy segment
				printf("copy segment %d %s\n",i,configinfile);

				#ifdef USE_MULTIGPU
					readasciivacconfig(configinfile,*p, meta, state, wmod,wd, hlines,mode);
				#else
					readbinvacconfig(configinfile,*p, meta, wmod,wd, *state );
				#endif
				gathersegment(*p, wnew,wdnew,wmod,wd);
				printf(" here read and gath nt=%d pid=%d i=%d\n",n,p->ipe,i);
			}

			p->n[0]=ni;
			p->n[1]=nj;
			#ifdef USE_SAC_3D
				p->n[2]=nk;
			#endif
			#ifdef USE_SAC_3D
				sprintf(configinfile,"%s%d.out",tcfg,n);  	     
			#else
				sprintf(configinfile,"%s%d.out",tcfg,n);  	     	     
			#endif           

			state->it=n;
			sprintf(configfile,"%s",cfggathout);
			writevacgatherconfig(configfile,n,*p, meta , wnew,wdnew,*state);
			printf(" here configfile %s nt=%d pid=%d \n",configfile,n,p->ipe);
			p->ipe=myipe;

		}//if p->ipe==0

		gpusync();
		//}//loop over nt steps
		//printf("proc %d here \n", p->ipe);

	}//if mode==gather
	#endif
	/*********************************************************************************************************/
	/* End of section to gather data
	/*********************************************************************************************************/


	/*********************************************************************************************************/
	/* Start of section to run special user initialisation
	/*********************************************************************************************************/
        //special user initialisation for the configuration 
        //this is a parallel routine
        if(mode==init)
        {
		p->mode=mode;

		#ifdef USE_MULTIGPU
			gpusync();
		#endif
                printf("init_config\n");
		initconfig(p, &meta, wmod,wd);
                printf("user initialisation\n");
		initialisation_user1(wmod,wd,p);

		// initialisation_user2(wmod,wd,p);
		//write the config file to ascii
                printf("writing ini file\n");
		writeasciivacconfig(configinfile,*p, meta , wmod,wd,hlines,*state,mode);
		#ifdef USE_MULTIGPU
			gpusync();
		#endif
        }   
	/*********************************************************************************************************/
	/* End of section to run special user initialisation
	/*********************************************************************************************************/




	//p->it=0;
	int order=0;
        
        if(mode==run)
        {
        //intialise arrays on GPU

	cuinit(&p,&bp,&wmod,&wnew,&wd,&state,&d_p,&d_bp,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
       /*********************************************************************************************************/
       /* Start of grid initialisation */
       /*********************************************************************************************************/

        //same as the grid initialisation routine in SAC
        //ensures boundaries defined correctly
	#ifdef USE_MPI

		cuinitmgpubuffers(&p, &w, &wmod, &temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,   &gmpiw0, &gmpiwmod0,   &gmpiw1, &gmpiwmod1,   &gmpiw2, &gmpiwmod2, &d_p, &d_w, &d_wmod,&d_wtemp2,  &d_gmpivisc0,  &d_gmpivisc1,  &d_gmpivisc2, &d_gmpiw0, &d_gmpiwmod0, &d_gmpiw1, &d_gmpiwmod1, &d_gmpiw2, &d_gmpiwmod2);


	/*cuinitmgpurbuffers(&p,    
			&d_gmpiviscr0,    
			&d_gmpiviscr1,    
			&d_gmpiviscr2,   
			&d_gmpiwr0, 
			&d_gmpiwmodr0,   
			&d_gmpiwr1, 
			&d_gmpiwmodr1,   
			&d_gmpiwr2, 
			&d_gmpiwmodr2);*/


		gpusync();
		int iii[3];
		int ip,jp;
		iii[2]=0;
		p->it=-1;
		
		printf("buffer initialisation complete %d\n",p->ipe);
gpusync();


		cucopywdtompiwd(&p,&wd,    &gmpiw0,     &gmpiw1,    &gmpiw2, &d_p,  &d_wd,    &d_gmpiw0,   &d_gmpiw1,   &d_gmpiw2,  order,0);
		gpusync();
				printf("%d here 1\n",p->ipe);

//still using hostcopy method temporarily
//#ifdef USE_GPUDIRECT
//		mpibound(NDERV, d_gmpiw0,d_gmpiw1,d_gmpiw2, p,0);
//		gpusync();
//		cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiwr0,    &d_gmpiwr1,    &d_gmpiwr2, order,0);
						printf("%d here 2\n",p->ipe);
//		gpusync();
//#else
		mpibound(NDERV, gmpiw0,gmpiw1,gmpiw2, p,0);
		gpusync();
		cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiw0,    &d_gmpiw1,    &d_gmpiw2, order,0);
		gpusync();
//#endif
				printf("%d here 3\n",p->ipe);
		cucopywdtompiwd(&p,&wd,    &gmpiw0,     &gmpiw1,    &gmpiw2, &d_p,  &d_wd,    &d_gmpiw0,   &d_gmpiw1,   &d_gmpiw2,  order,1);
		gpusync();
//#ifdef USE_GPUDIRECT
//				printf("%d here 4\n",p->ipe);
//		mpibound(NDERV, d_gmpiw0,d_gmpiw1,d_gmpiw2 , p,1);
//		gpusync();
//		cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiwr0,    &d_gmpiwr1,    &d_gmpiwr2, order,1);
						printf("%d here 5\n",p->ipe);
//#else

                printf("call mpibound %d\n",p->ipe);
		mpibound(NDERV, gmpiw0,gmpiw1,gmpiw2 ,p,1);
	        printf("leave mpibound %d\n",p->ipe);
		gpusync();
		cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiw0,    &d_gmpiw1,    &d_gmpiw2, order,1);
//#endif

				printf("%d here 6\n",p->ipe);


#ifdef USE_SAC_3D
			gpusync();
			cucopywdtompiwd(&p,&wd,    &gmpiw0,     &gmpiw1,    &gmpiw2, &d_p,  &d_wd,    &d_gmpiw0,   &d_gmpiw1,   &d_gmpiw2,  order,2);
			gpusync();
	#ifdef USE_GPUDIRECT
			mpibound(NDERV, d_gmpiw0,d_gmpiw1,d_gmpiw2  , p,2);
			gpusync();
			cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   d_gmpiwr0,    d_gmpiwr1,    d_gmpiwr2, order,2);

	#else
			mpibound(NDERV, gmpiw0,gmpiw1,gmpiw2, p,2);
			gpusync();
			cucopywdfrommpiwd(&p,&wd,     &gmpiw0,     &gmpiw1,     &gmpiw2,  &d_p,  &d_wd,   &d_gmpiw0,    &d_gmpiw1,    &d_gmpiw2, order,2);
	#endif
#endif






		p->it=n+1;
		printf("1 after buffer copy %d\n",p->ipe);
		cuupdatehostwd(&p,&wd,&wmod,&temp2,&state,&d_p,&d_wd,&d_wmod,&d_wtemp2,  &d_state,n);
		printf("2 after buffer copy %d\n",p->ipe);
		initgrid(&p,&state,&wd,&d_p, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
		printf("grid initialised\n");
	        cusync(&p);	   
	  #endif     //endif MPI

	cusync(&p);
         #ifdef USE_MPI
		p->it=n+1;
		gpusync();
		cusync(&p);
        #endif //use_MPI


      #ifndef USE_MULTIGPU
	int iii[2];
	initgrid(&p,&state,&wd,&d_p, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
      #endif
       /*********************************************************************************************************/
       /* End of grid initialisation */
       /*********************************************************************************************************/

	//mgpuinit_stage2(p);

       /*********************************************************************************************************/
       /* Start of section initialising boundaries */
       /*********************************************************************************************************/
      
	for(int ii=0; ii<NVAR; ii++)
	for(int idir=0; idir<NDIM; idir++)
        for(int ibound=0; ibound<2; ibound++)
	{
	   p->it=-1;  //initialise fixed boundaries
	   if((p->boundtype[ii][idir][ibound])==5)  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
	   {
		       cuboundary(&p, &bp, &d_p, &d_bp, &d_state, &d_wmod, 0,idir,ii);	 
	   }
	}
       /*********************************************************************************************************/
       /* End of section initialising boundaries */
       /*********************************************************************************************************/
       
        p->it=its;  
        



       /*********************************************************************************************************/
       /* Start of section initialising the configuration */
       /*********************************************************************************************************/

        #ifdef USE_MPI

        for(int ordert=0; ordert<=1; ordert++)//ozt case may only work with case ordert==1
        for(int idir=0; idir<NDIM;idir++)
        {
		//for runge kutta will need to run this several times  for each order 
		if(p->ipe==0)          
		printf("before mpi trans mpiwmod\n");
          // if(idir==1)
		 cucopywtompiwmod(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, ordert,idir);

		gpusync();
		if(p->ipe==0)          
		printf("mpi trans mpiwmod\n");
		
#ifdef USE_GPUDIRECT
		mpiboundmod(NVAR, d_gmpiwmod0,d_gmpiwmod1,d_gmpiwmod2 ,p,idir);
#else
		mpiboundmod(NVAR, gmpiwmod0,gmpiwmod1,gmpiwmod2 ,p,idir);
#endif
		gpusync();
		//for runge kutta will need to run this several times  for each order  
           //if(idir==1) 
		if(p->ipe==0)          
		printf("after mpi trans mpiwmod\n");
#ifdef USE_GPUDIRECT
		cucopywmodfrommpiw(&p,&w, &wmod,      &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmodr0,   &d_gmpiw1, &d_gmpiwmodr1,   &d_gmpiw2, &d_gmpiwmodr2,ordert,idir);
#else      
		cucopywmodfrommpiw(&p,&w, &wmod,      &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,ordert,idir);
#endif


		gpusync();		



         }




	#endif






       /*********************************************************************************************************/
       /* End of section initialising the configuration */
       /*********************************************************************************************************/




       /*********************************************************************************************************/
       /* Start of section getting parameters from metadata file */
       /*********************************************************************************************************/

	//For a steerable simulation generate and save a dxformfile that saves a single data step
	//used for the steering dx module
	//printf("here in runsim2a\n");
	#ifdef USE_IOME
	getmetadata_(elist.id,"directory",&sdir,elist.port,elist.server);
	//sdir=metadata.directory
	//name=metadata.name;
	getmetadata_(elist.id,"name",&name,elist.port,elist.server);
	//disp(sdir,name)
	//printf("here in runsim3\n");
	sprintf(outfile,"%s/%s.out",sdir,name);
	#endif
        /*********************************************************************************************************/
        /* End of section getting parameters from metadata file */
        /*********************************************************************************************************/







        /*********************************************************************************************************/
        /* Start of section to set steering control (n.b. commented out) */
        /*********************************************************************************************************/     
	//createlog(meta.log_file);
	//while(finishsteering == 0)
	//{
	 
	  //  if( steeringenabled==0)
	  //    finishsteering=1;
        /*********************************************************************************************************/
        /* End of section to enable steering control */
        /*********************************************************************************************************/     







	
	real t1,t2,ttot;
	int ordero=0;
	int order1;
	int orderb=0;
	int ii,ii0,ii1;
	real dtdiffvisc,dtgrav,dttemp,ga;
	ttot=0;
	real time=0.0;
	state->it=0;
	state->t=0;
	state->dt=p->dt;



       /*********************************************************************************************************/
       /* Apply boundaries */
       /*********************************************************************************************************/
	for(int ii=0; ii<=(b1+(NDIM-1)); ii++)
	for(int idir=0; idir<NDIM; idir++)
        //for(int ibound=0; ibound<2; ibound++)
	{
	   

		       cuboundary(&p, &bp, &d_p, &d_bp, &d_state, &d_wmod, 1,idir,ii);	 

	}


       /*********************************************************************************************************/
       /* Start looping over iterations*/
       /*********************************************************************************************************/
        printf("its %d\n",nt);
        //for(n=nt+1;n<=nt;n++)
	for( n=its;n<=nt;n++)
	{
	    	p->it=n;	
		//gpusync();
		//printf("%d,step %d\n",p->ipe,n);
		//gpusync();
	
		if((p->rkon)==0)
		{
	  		cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,0);
        	}	
	

	    if(((n-1)%(p->cfgsavefrequency))==0)
	    {
			//writeconfig(name,n,*p, meta , w);

                //note below writing wmod to output no need for the w field
                //note the zeroth component of wmod written 
		#ifndef USE_MPI
			// writevtkconfig(configfile,n,*p, meta , w);
                     // writevacconfig(configfile,n,*p, meta , w,wd,*state);
                      writevacconfig(configfile,n,*p, meta , wmod,wd,*state);
		#else
		   //  writeasciivacconfig(configfile,*p, meta , w,wd,hlines,*state,mode);
                    writeasciivacconfig(configfile,*p, meta , wmod,wd,hlines,*state,mode);
                  ;//  writevacconfig(configfile,n,*p, meta , wmod,wd,*state);
                 #endif
		//writevacconfig(configfile,n,*p, meta , wmod,wd,*state);


	 
                printf("finished write routine\n");
	    }




	

	    order=0;
	    t1=second();	
	





       /*********************************************************************************************************/
       /* Start single step  iteration*/
       /*********************************************************************************************************/
	   if(p->moddton==1.0)
	   {


                tc=second();
		p->maxcourant=0.0;
                dtgrav=BIGDOUBLE;
		courantmax=SMALLDOUBLE;
              dt=0.005;
              p->dt=0.005;
		for(int dim=0; dim<=(NDIM-1); dim++)
		{
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			//cucomputemaxcourant(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);  //potential bottleneck here
		}

		for(int dim=0; dim<=(NDIM-1); dim++)
		{			
                           dttemp=(p->cmax/(p->dx[dim]));
                           if(dttemp>courantmax ) courantmax=dttemp;
                }
               p->maxcourant=courantmax;


		#ifdef USE_MPI
                   tv=second();
                   gpusync();
		   mpiallreduce(&(p->maxcourant), MPI_MAX);
                   tcom+=(second()-tv);
		#endif
	
		if(     (dttemp=(  (p->courant)/(p->maxcourant)  ))>SMALLDOUBLE  && dttemp<dt  )
		       p->dt=dttemp;
		printf("new dt is %g %g %g\n",(p->courant)/(p->maxcourant),p->dt,(p->maxcourant));

		//if(n>1)
		//   cugetdtvisc1(&p,&d_p,&d_wmod, &wd,&d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);


		//include hyperdiffusion contribution
                if(n>(its+1))
                {
                dtdiffvisc=BIGDOUBLE;
		for(int dim=0; dim<=(NDIM-1); dim++)
		{			
                           dttemp=0.25/(p->maxviscoef/(p->dx[dim]));
                           if(dttemp<dtdiffvisc && dttemp>SMALLDOUBLE) dtdiffvisc=dttemp;
                }

                ;//if(dtdiffvisc<(p->dt) && dtdiffvisc>SMALLDOUBLE) p->dt=dtdiffvisc;
                }


                tcal+=(second()-tc);
               printf(" dtdiffvisc %20.10g %20.10g  %20.10g\n",dttemp,p->maxviscoef,p->dtdiffvisc);

                //printf("ipe %d dtdiffvisc %20.10g  %20.10g\n",p->ipe,p->maxviscoef,p->dtdiffvisc);
		#ifdef USE_MPI
                   tv=second();
                   gpusync();
		   mpiallreduce(&(p->dtdiffvisc), MPI_MIN);
                   tcom+=(second()-tv);
		#endif

		
		//if((p->dtdiffvisc)>SMALLDOUBLE && (p->dt)>((p->dtdiffvisc)) )
		//	                      			p->dt=(p->dtdiffvisc);
		#ifdef USE_MPI
			printf(" on pe %d modified dt is %20.10g \n",p->ipe,p->dt);
		#else
			printf("modified dt is %20.10g  %20.10g\n",p->dt,p->dtdiffvisc);
		#endif
		//include gravitational modification
		/*for(int dim=0; dim<=(NDIM-1); dim++)
		{			
			if((ga=abs(p->g[dim])) > 0)
                        {
                           dttemp=1.0/sqrt(ga/(p->dx[dim]));
                           if(dttemp<dtgrav && dttemp>SMALLDOUBLE) dtgrav=dttemp;
                         }
                }

                if(dtgrav<(p->dt) && dtgrav>SMALLDOUBLE) p->dt=dtgrav;*/
                p->maxviscoef=SMALLDOUBLE;




	   } 
       /*********************************************************************************************************/
       /* End of single step  iteration*/
       /*********************************************************************************************************/








       /*********************************************************************************************************/
       /* Start single step  iteration*/
       /*********************************************************************************************************/
	if((p->rkon)==0)
	{
	  ordero=1;
	  order=0;
         tc=second();
         p->qt=(p->qt)+dt;
	 // cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
	  
	 for(int dir=0;dir<NDIM; dir++)
	 {
		 cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir); 
		      
         }
	 for(int dir=0;dir<NDIM; dir++)
	 {
		  
		 for(int f=rho; f<=(mom1+NDIM-1); f++)
		  { 

                  #ifdef USE_SAC_3D
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  || (f==mom3 && dir==2) )
                  #else
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  )
                  #endif
                    {
		       cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
                       cucomputepbg(&p,&d_p,&d_wmod, &d_wd,order,dir);
                     }
		     cucentdiff1(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);	     
		  } //end looping over fields for cucentdiff1
		   #ifndef ADIABHYDRO
		   for(int f=energy; f<=(b1+(NDIM-1)); f++)
		   {
		     if(f==energy)
		     {
			 cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
			 cucomputepbg(&p,&d_p,&d_wmod, &d_wd,ordero,dir);
			 cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
		     }	      
		     cucentdiff2(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
		   }//end looping over fields for cucentdiff2
		   #endif
	  }//end loop over directions
	   
	 
	  cugrav(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);//gravitational contributions

	  if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt);
          tcal+=(second()-tc);

           /*********************************************************************************************************/
           /* Start  of hyperdiffusion contributions for single step*/
           /*********************************************************************************************************/
	   if(p->hyperdifmom==1)
	   {
            tc=second();
	    dt=(p->dt);
		     p->maxviscoef=0.0;
	    
	    #ifdef USE_SHOCKVISC
	       cunushk1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
	    #endif
            tcal+=(second()-tc);
            //density hyperdiffusion term
	    for(int dim=0; dim<=(NDIM-1); dim++)
	    {
			tc=second();		      
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                     //p->cmax=2.0;
			 tcal+=(second()-tc);
		      #ifdef USE_MPI
			      tv=second();
                              gpusync();
			      mpiallreduce(&(p->cmax), MPI_MAX);

                              tcom+=(second()-tv);
		      #endif
		      cmax[dim]=p->cmax;
		      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

		      #ifdef USE_MPI

                          tv=second();
			  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);

                          gpusync();
			  mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
                          gpusync();

			  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			  tcom+=(second()-tv);

		      #endif
                      tc=second();
		      cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
		      cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

	              cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,p->dt);
                      tcal+=(second()-tc);
	     } //end of rho hyperdif contributions for each direction


             //energy hyperdiffusion term
	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
	        //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
		//cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                tc=second();
		p->cmax=cmax[dim];
		#ifdef USE_MPI
		      ;//mpiallreduce(&(p->cmax), MPI_MAX);
		#endif
	        cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
                tcal+=(second()-tc);
		#ifdef USE_MPI
                        tv=second();
			cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
			cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                        tcom+=(second()-tv);
		#endif
                tc=second();
		cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
		cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);	       
	        cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt); 
                tcal+=(second()-tc);  
	     }

        //momentum hyperdiffusion term
	for(int dim=0; dim<=(NDIM-1); dim++)
	       for(int f=0; f<=(NDIM-1); f++)		   	                 
		{
			  //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			  //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                          p->cmax=cmax[dim];
			  #ifdef USE_MPI
			     ;// mpiallreduce(&(p->cmax), MPI_MAX);
			  #endif
                          tc=second();
		          cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                          tcal+=(second()-tc);
		          #ifdef USE_MPI
				  tv=second();
                                  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				 mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                        tcom+=(second()-tv);
	                 #endif
                          tc=second();
			 cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			 cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
tcal+=(second()-tc);

tc=second();
		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=f;
		                  }
		                  else
		                  {
				           ii=f;
				           ii0=dim;
		                   }

				  if(ii==dim)
				    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
				  else
				    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
		        }
                        tcal+=(second()-tc);
		     }  //end of loop over dimensions and fields


                     //b field hyperdiffusion term
		     int jj,mm,kk;
		     real sb;
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++) 
		     if(f!=dim)           
		     {
			       //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			       //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			       p->cmax=cmax[dim];
			       #ifdef USE_MPI
			      		;//mpiallreduce(&(p->cmax), MPI_MAX);
			       #endif
                          tc=second();
			       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
                tcal+=(second()-tc);
			       #ifdef USE_MPI
                                  tv=second();
				  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				  mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);	
                        tcom+=(second()-tv);	 
		               #endif
                          tc=second();
			       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			       cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);

                                
				for(ii1=0;ii1<=1;ii1++)  //start of compute cross-product term of b field hyperdiffusion terms
				{
				          if (ii1 == 0)
				          {
						   jj=dim;
						   mm=f;
						   sb=-1.0;
						   ii0=dim;
				          }
				          else
				          {
						   ii0=f;
						   mm=dim;
						   sb=1.0;
						   jj=f;				           
				          }

					  if(mm==dim)
					     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
					  else
					     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);		
				}//end of compute cross-product term of b field hyperdiffusion terms
                              tcal+=(second()-tc);
		      }//end of loop over fields and dimensions

	   }//closes if(p->hyperdifmom==1)
           /*********************************************************************************************************/
           /* End of hyperdiffusion contributions for single step */
           /*********************************************************************************************************/

	  //source terms
          tc=second();
          cusource(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);
                tcal+=(second()-tc);

          //HFFILTER TERM
          /*if(((n-1)%(p->hffiltfrequency))==0)
          {
          	tc=second();
          	cuhffilt(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);
          	tcal+=(second()-tc);
          }*/

	for(int ii=0; ii<=(b1+(NDIM-1)); ii++)
	for(int idir=0; idir<NDIM; idir++)
        //for(int ibound=0; ibound<2; ibound++)
	  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,idir,ii);

	} //end of if((p->rkon)==0)
       /*********************************************************************************************************/
       /* End single step  iteration*/
       /*********************************************************************************************************/






       /*********************************************************************************************************/
       /* Start runge-kutta  iteration*/
       /*********************************************************************************************************/

	   if((p->rkon)==1) p->maxviscoef=0.0;

	   if((p->rkon)==1)
	   for(order=0; order<4; order++) 
	   {	   
		   ordero=order+1;
		   dt=(p->dt)/2.0;
		   orderb=order+2;

		   if(order==2)
		   {
		      dt=(p->dt);
		      orderb=1;
		    }


		   if(order==3)
		   {
		      dt=(p->dt)/6.0;
		      ordero=0;
		      orderb=0;
		   }


		   //cucomputedervfields(&p,&d_p,&d_wmod, &d_wd,order);
                   for(int dir=0;dir<(NDIM-1); dir++)
		   {
			cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
                   }




	      for(int dir=0;dir<(NDIM-1); dir++)
	     {
		 for(int f=rho; f<=(mom1+NDIM-1); f++)
		  { 

                  #ifdef USE_SAC_3D
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  || (f==mom3 && dir==2) )
                  #else
		      if((f==mom1 && dir==0)  ||  (f==mom2 && dir==1)  )
                  #endif
                    {
		       cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
                       cucomputepbg(&p,&d_p,&d_wmod, &d_wd,order,dir);
                     }
		     cucentdiff1(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order,ordero,p->dt,f,dir);	     
		  } //end looping over fields for cucentdiff1


		   #ifndef ADIABHYDRO
		   for(int f=energy; f<=(b1+(NDIM-1)); f++)
		   {
		     if(f==energy)
		     {
			 cucomputevels(&p,&d_p,&d_wmod, &d_wd,order,dir);
			 cucomputepbg(&p,&d_p,&d_wmod, &d_wd,ordero,dir);
			 cucomputept(&p,&d_p,&d_wmod, &d_wd,order,dir);
		     }	      
		     cucentdiff2(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt,f,dir);
		   }//end looping over fields for cucentdiff2
		   #endif
           }//end of loop over dimensions



		   if(p->divbon==1)
		       cudivb(&p,&d_p,&d_w,&d_wmod, &d_dwn1, &d_wd, order,ordero,p->dt);


	           cugrav(&p,&d_p,&d_state,&d_wmod,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);//gravitational contributions



		   /*********************************************************************************************************/
		   /* End of hyperdiffusion contributions for multi step */
		   /*********************************************************************************************************/
		   if(p->hyperdifmom==1)
		   {



            tc=second();
	    dt=(p->dt);

	    
	    #ifdef USE_SHOCKVISC
	       cunushk1(&p,&d_p,&d_wmod, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2);
	    #endif
            tcal+=(second()-tc);
            //density hyperdiffusion term
	    for(int dim=0; dim<=(NDIM-1); dim++)
	    {
			tc=second();		      
			cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                     //p->cmax=2.0;
			 tcal+=(second()-tc);
		      #ifdef USE_MPI
			      tv=second();
                              gpusync();
			      mpiallreduce(&(p->cmax), MPI_MAX);

                              tcom+=(second()-tv);
		      #endif
		      cmax[dim]=p->cmax;
		      cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

		      #ifdef USE_MPI

                          tv=second();
			  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);

                          gpusync();
			  mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
                          gpusync();

			  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			  tcom+=(second()-tv);

		      #endif
                      tc=second();
		      cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);
		      cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,rho,dim);

	              cuhyperdifrhosource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,rho,dim,p->dt);
                      tcal+=(second()-tc);
	     } //end of rho hyperdif contributions for each direction


             //energy hyperdiffusion term
	     for(int dim=0; dim<=(NDIM-1); dim++)
	     {
	        //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
		//cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                tc=second();
		p->cmax=cmax[dim];
		#ifdef USE_MPI
		      ;//mpiallreduce(&(p->cmax), MPI_MAX);
		#endif
	        cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
                tcal+=(second()-tc);
		#ifdef USE_MPI
                        tv=second();
			cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
			mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
			cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                        tcom+=(second()-tv);
		#endif
                tc=second();
		cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);
		cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,energy,dim);	       
	        cuhyperdifesource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,energy,dim,dt); 
                tcal+=(second()-tc);  
	     }

        //momentum hyperdiffusion term
	for(int dim=0; dim<=(NDIM-1); dim++)
	       for(int f=0; f<=(NDIM-1); f++)		   	                 
		{
			  //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			  //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
                          p->cmax=cmax[dim];
			  #ifdef USE_MPI
			     ;// mpiallreduce(&(p->cmax), MPI_MAX);
			  #endif
                          tc=second();
		          cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
                          tcal+=(second()-tc);
		          #ifdef USE_MPI
				  tv=second();
                                  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				 mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
                        tcom+=(second()-tv);
	                 #endif
                          tc=second();
			 cuhyperdifvisc1r(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
			 cuhyperdifvisc1l(&p,&d_p,&d_wmod,&wd,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,mom1+f,dim);
tcal+=(second()-tc);

tc=second();
		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=f;
		                  }
		                  else
		                  {
				           ii=f;
				           ii0=dim;
		                   }

				  if(ii==dim)
				    cuhyperdifmomsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
				  else
				    cuhyperdifmomsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,ii,ii0,p->dt);
		        }
                        tcal+=(second()-tc);
		     }  //end of loop over dimensions and fields


                     //b field hyperdiffusion term
		     int jj,mm,kk;
		     real sb;
		     for(int dim=0; dim<=(NDIM-1); dim++)
		     for(int f=0; f<=(NDIM-1); f++) 
		     if(f!=dim)           
		     {
			       //cucomputec(&p,&d_p,&d_wmod, &d_wd,order,dim);
			       //cucomputemaxc(&p,&d_p,&d_wmod, &d_wd,order,dim,&wd,&d_wtemp);
			       p->cmax=cmax[dim];
			       #ifdef USE_MPI
			      		;//mpiallreduce(&(p->cmax), MPI_MAX);
			       #endif
                          tc=second();
			       cuhyperdifvisc1ir(&p,&d_p,&d_wmod,  &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
                tcal+=(second()-tc);
			       #ifdef USE_MPI
                                  tv=second();
				  cucopytompivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);
				  mpivisc(dim,p,gmpivisc0,gmpivisc1,gmpivisc2);
				  cucopyfrommpivisc(&p,&temp2, &gmpivisc0, &gmpivisc1, &gmpivisc2,  &d_p,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2);	
                        tcom+=(second()-tv);	 
		               #endif
                          tc=second();
			       cuhyperdifvisc1r(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);
			       cuhyperdifvisc1l(&p,&d_p,&d_wmod, &wd, &d_wd,order,&d_wtemp,&d_wtemp1,&d_wtemp2,b1+f,dim);

                                
				for(ii1=0;ii1<=1;ii1++)  //start of compute cross-product term of b field hyperdiffusion terms
				{
				          if (ii1 == 0)
				          {
						   jj=dim;
						   mm=f;
						   sb=-1.0;
						   ii0=dim;
				          }
				          else
				          {
						   ii0=f;
						   mm=dim;
						   sb=1.0;
						   jj=f;				           
				          }

					  if(mm==dim)
					     cuhyperdifbsource1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);
					  else
					     cuhyperdifbsourcene1(&p,&d_p,&d_wmod, &d_dwn1, &d_wd,order,ordero,&d_wtemp,f,dim,jj,ii0,mm,sb,p->dt);		
				}//end of compute cross-product term of b field hyperdiffusion terms
                              tcal+=(second()-tc);
		      }//end of loop over fields and dimensions





		   }//closes if(p->hyperdifmom==1)



	  //source terms
          tc=second();
          cusource(&p,&d_p,&d_state,&d_w,&d_wmod, &d_dwn1, &d_wd,order, ordero,p->dt);
                tcal+=(second()-tc);



		   cuadvance(&p,&d_p,&d_wmod,&d_w,order);
		   
		   
		   
		   #ifdef USE_MPI
                           for(int idir=0; idir<NDIM;idir++)
        {
			cucopywtompiwmod(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order,idir);


		  // mpibound(NVAR, gmpiw0,gmpiw1,gmpiw2 ,p,idir);
		#ifdef USE_GPUDIRECT
		   mpiboundmod(NVAR, d_gmpiwmod0,d_gmpiwmod1,d_gmpiwmod2,p,idir);
		#else
		   mpiboundmod(NVAR, gmpiwmod0,gmpiwmod1,gmpiwmod2,p,idir);
		#endif


		#ifdef USE_GPUDIRECT

			cucopywmodfrommpiw(&p,&w, &wmod,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmodr0,   &d_gmpiw1, &d_gmpiwmodr1,   &d_gmpiw2, &d_gmpiwmodr2,order,idir);

		#else


			cucopywmodfrommpiw(&p,&w, &wmod,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,    &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order,idir);
		#endif
      }

		   #endif
		   
	for(int ii=0; ii<=(b1+(NDIM-1)); ii++)
	for(int idir=0; idir<NDIM; idir++)
        //for(int ibound=0; ibound<2; ibound++)
	  cuboundary(&p,&bp,&d_p,&d_bp,&d_state,&d_wmod, ordero,idir,ii);		   
		   

		   

	   }//looping over orders
       /*********************************************************************************************************/
       /* End runge-kutta  iteration*/
       /*********************************************************************************************************/


	  #ifdef USE_MPI




   // if((p)->ipe==0 && ((p)->it)==2)
   // {
  //       printf("ipe3 mpiwmod \n");
    //     for(int i=0; i<4*((p)->n[0]);i++)
    //         printf("%d %lg \n",i, (gmpiwmod1[i]));
         //printf("\n");


        /* printf("ipe3 mpiwmod \n");
         for(int i=0; i<4*((p)->n[0]);i++)
             printf("%d %lg ",i,(gmpiw0[i]));
         printf("\n");*/
    // }
        order=1;
       
        tcom1=second();
        for(int idir=0; idir<NDIM;idir++)
        {




                  gpusync();
		 //  cucopywtompiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order,idir);
              //if(idir==1)


		   cucopywtompiwmod(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2, order,idir);

                 gpusync();
	      //   mpibound(NVAR, gmpiw0,gmpiw1,gmpiw2 ,p,idir);
		 //  gpusync();
                #ifdef USE_GPUDIRECT
		   mpiboundmod(NVAR, d_gmpiwmod0,d_gmpiwmod1,d_gmpiwmod2 ,p,idir);

                #else

		   mpiboundmod(NVAR, gmpiwmod0,gmpiwmod1,gmpiwmod2, p,idir);
                #endif

		 
 gpusync();

		//   cucopywfrommpiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order,idir);	
		// if(idir==1)
                #ifdef USE_GPUDIRECT
		   cucopywmodfrommpiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmodr0,   &d_gmpiw1, &d_gmpiwmodr1,   &d_gmpiw2, &d_gmpiwmodr2,order,idir);	
                #else



		   cucopywmodfrommpiw(&p,&w, &wmod,    &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,  &d_w, &d_wmod,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2,order,idir);	
		#endif


 gpusync();






   }
     tcom2=second()-tcom1;
     tcom+=tcom2;
	   
	  #endif





	   p->it=n+1;

        //initgrid(&p,&w,&wnew,&state,&wd,&d_p,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
        //cuupdate(&p,&w,&wmod,&temp2,&state,&d_gp[igid],&d_gw[igid],&d_gwmod[igid],&d_gwtemp2[igid],  &d_gstate[igid],n);
        //if(igid != 3)
         tc=second();
         cuupdate(&p,&w,&wmod,&temp2,&state,&d_p,&d_w,&d_wmod,&d_wtemp2,  &d_state,n);
                tcal+=(second()-tc);
	//initgrid(&p,&w,&wnew,&state,&wd,&d_p,&d_gw[igid],&d_wnew,&d_wmod, &d_dwn1,  &d_gwd[igid], &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
	//initgrid(&p,&w,&wnew,&state,&wd,&d_gp[igid],&d_gw[igid],&d_gwnew[igid],&d_gwmod[igid], &d_gdwn1[igid],  &d_gwd[igid], &d_gstate[igid],&d_gwtemp[igid],&d_gwtemp1[igid],&d_gwtemp2[igid]);

       // igid=0;
        cusync(&p);











	   t2=second()-t1;
	   ttot+=t2;
	   printf("step %d time total=%f com=%f cal=%f\n",n,ttot,tcom,tcal);

	   state->it=n;
	   state->t=time+(p->dt);
	   time=state->t;
	   state->dt=p->dt;



	   //appendlog(meta.log_file,*p, *state);
	    /*getintparam_(&elist.id,"steeringenabled",&steeringenabled,&elist.port,elist.server);
	    if(steeringenabled==1)
	    {
	      //disp('getting updatea params');
	      //for steering get the modified control params
	      double dg;
	      getintparam_(&elist.id,"finishsteering",&finishsteering,&elist.port,elist.server);//source y location  
		// Constants
	      getdoubleparam_(elist.id,"g",&dg,elist.port,elist.server);

	      g=dg;
	     
	    }*/

           // gpusync();





	
	
	
	}
       /*********************************************************************************************************/
       /* End of looping over iterations*/
       /*********************************************************************************************************/

         //cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
         //printf("at cufinish end here %d\n",p->ipe);
         printf("at cufinish end here\n");


        

        //igid=0;
        cusync(&p);


	
        } //mode=0 clean up routine//if(mode==run)









//}//if(mode==run)

cusync(&p);

//cufinish(&p,&w,&wnew,&state,&d_p,&d_bp,&d_w,&d_wnew,&d_wmod, &d_dwn1,  &d_wd, &d_state,&d_wtemp,&d_wtemp1,&d_wtemp2);
#ifdef USE_MPI
			gpusync();
	     printf("at cumpifinish end here %d\n",p->ipe);
#endif
	free(hlines);
	free(p);
	free(bp);
	free(sdir);
	free(name);
	free(outfile);
	free(formfile);

	#ifdef USE_MPI
				gpusync();
	    // printf("at cumpifinish end here %d\n",p->ipe);
             mgpufinalize(p);

	   ;// cufinishmgpu(&p,&w, &wmod, &temp2,&gmpivisc0,&gmpivisc1,&gmpivisc2,   &gmpiw0, &gmpiwmod0,    &gmpiw1, &gmpiwmod1,    &gmpiw2, &gmpiwmod2, &d_p,   &d_w, &d_wmod,&d_wtemp2,    &d_gmpivisc0,    &d_gmpivisc1,    &d_gmpivisc2,   &d_gmpiw0, &d_gmpiwmod0,   &d_gmpiw1, &d_gmpiwmod1,   &d_gmpiw2, &d_gmpiwmod2);
           ;// MPI::Finalize();

	#endif
            printf("return\n");
		return 0;
	

}


int fencode3_test (struct params *dp,int *ii, int field) {

#ifdef USE_SAC_3D
   return (ii[2]*((dp)->n[0])*((dp)->n[1])  + ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
#else
   return ( ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])));
#endif

}


int encode3_in(struct params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

void initconfig(struct params *k, struct meta *md, real *w, real *wd)
{




	int i1,j1,k1;
        int ni=k->n[0];
        int nj=k->n[1];
        unsigned long int ilv;
        printf("%d %d\n",ni,nj);
        for(i1=0; i1<(k->n[0]) ;i1++)
	  for(j1=0; j1<(k->n[1]) ;j1++)
          {
                    
                    wd[encode3_in(k,i1,j1,k1,pos1)]=(k->xmin[0])+((real)i1)*((k->xmax[0])- (k->xmin[0])  )/ni;
                    wd[encode3_in(k,i1,j1,k1,delx1)]=((k->xmax[0])- (k->xmin[0])  )/ni;
                    wd[encode3_in(k,i1,j1,k1,pos2)]=(k->xmin[1])+((real)j1)*((k->xmax[1])- (k->xmin[1])  )/nj;
                    wd[encode3_in(k,i1,j1,k1,delx2)]=((k->xmax[1])- (k->xmin[1])  )/nj;

		    #ifdef USE_SAC3D
		            wd[encode3_in(k,i1,j1,k1,pos3)]=(k->xmin[2])+((real)k1)*((k->xmax[2])-(k->xmin[2]))/nk;
		            wd[encode3_in(k,i1,j1,k1,delx3)]=((k->xmax[2])- (k->xmin[2])  )/nk;
                    #endif

                    for(int f=rho; f<NVAR; f++)
                    //for(int f=rho; f<=b2; f++)
                    {
                    ;//ilv=j1*ni+i1+(ni*nj*f);
                     ilv=j1*ni+i1+(ni*nj*f);
                    w[ilv]=0.0;
                    switch(f)
		            {
		              case rho:
		            	w[ilv]=1.0;
			      break;
		              case mom1:
		            	w[ilv]=0.01;
			      break;
		              case mom2:
		            	w[ilv]=0.01;
			      break;
		              //case mom3:
		            	//w[j1*ni+i1+(ni*nj*f)]=0.0;
			      //break;
		              case energy:
		            	w[ilv]=0.0;
			      break;
		              case b1:
		            	w[ilv]=0.0;
			      break;
		              case b2:
		            	w[ilv]=0.0;
			      break;
		              //case b3:
		            //	w[j1*ni+i1+(ni*nj*f)]=0.0;
			     // break;
		            }; //end of switch to check for field

			}//end of loop over f

          }//end of loop over j and i
          printf("ilv %ld\n", ilv);


}



#ifdef USE_MPI
/*        #include "smaugmpi.h"*/


//#include "mpi.h"
//#include "iotypes.h"

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>





int sacencodempivisc0 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (  bound* ((   ((p->n[1])+2)*((p->n[2])+2)   ))+       (    (     iy+iz*((p->n[1])+2)    )    )      );
  #else
    return (   bound*(    ((p->n[1])+2)  )  +   iy     );
  #endif
}



int sacencodempivisc1 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (bound* ((   ((p->n[0])+2)*((p->n[2])+2))      )+   ((ix+iz*((p->n[0])+2)))  );
  #else
    return (   bound*(    ((p->n[0])+2)  )  +   ix     );
  #endif

  return 0;
}


int sacencodempivisc2 (struct params *p,int ix, int iy, int iz, int bound,int dim) {
  #ifdef USE_SAC_3D
    return (bound* ((   ((p->n[0])+2)*((p->n[1])+2))      )+   (  (ix+iy*((p->n[0])+2))    ));
  #endif
  return 0;
}

int sacencodempiw0 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[1])*(p->n[2]))   )+bound*(            +  ((p->n[1])*(p->n[2]))      )+   (  (iy+iz*(p->n[1]))    ));
  #else
    return (4*field*(p->n[1]) +bound*((p->n[1]))  +   (iy));
  #endif
  return 0;
}

int sacencodempiw1 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[0])*(p->n[2]))   )+
bound*(            +  ((p->n[0])*(p->n[2]))      )+   (  (ix+iz*(p->n[0]))    ));
  #else
    return (4*field*(p->n[0]) +bound*((p->n[0]))  +   (ix));
  #endif
return 0;
}


int sacencodempiw2 (struct params *p,int ix, int iy, int iz, int field,int bound) {
  #ifdef USE_SAC_3D
    return (4*field*(         ((p->n[0])*(p->n[1]))   )+
bound*(            +  ((p->n[0])*(p->n[1]))      )+   (  (ix+iy*(p->n[0]))    ));
  #endif
   return 0;
}

int encode3p2_sacmpi (struct params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*(((dp)->n[0])+2)*(((dp)->n[1])+2)  + iy * (((dp)->n[0])+2) + ix)+(field*(((dp)->n[0])+2)*(((dp)->n[1])+2)*(((dp)->n[2])+2)));
  #else
    return ( (iy * (((dp)->n[0])+2) + ix)+(field*(((dp)->n[0])+2)*(((dp)->n[1])+2)));
  #endif
  return 0;
}

/*void mpiinit(params *p);
void mpifinalize(params *p);
void mpisetnpediped(params *p, char *string);
void ipe2iped(params *p);
void iped2ipe(params *p);*/

//!=============================================================================
//subroutine mpiinit
//
//! Initialize MPI variables
//

//!----------------------------------------------------------------------------
//call MPI_INIT(ierrmpi)
//call MPI_COMM_RANK (MPI_COMM_WORLD, ipe, ierrmpi)
//call MPI_COMM_SIZE (MPI_COMM_WORLD, npe, ierrmpi)

//! unset values for directional processor numbers
//npe1=-1;npe2=-1;
//! default value for test processor
//ipetest=0
void mgpuinit(struct params *p)
{
    int nmpibuffer;
   int numbuffers=2;
   int i;


     //MPI::Intracomm comm;
     //MPI_Init(&argc, &argv);

     /*gwall_time = MPI_Wtime();
     comm=MPI::COMM_WORLD;
     p->npe=comm.Get_size();
     p->ipe=comm.Get_rank();*/

#ifdef USE_SAC3D
   if((p->n[0])>=(p->n[1])  && (p->n[0])>=(p->n[2]))
   {
     if((p->n[1])>(p->n[2]))
       nmpibuffer=NDERV*(p->n[0])*(p->n[1])*(p->ng[0]);
     else
       nmpibuffer=NDERV*(p->n[0])*(p->n[2])*(p->ng[0]);
   }
   else if((p->n[1])>=(p->n[0])  && (p->n[1])>=(p->n[2]))
   {
     if((p->n[0])>(p->n[2]))
       nmpibuffer=NDERV*(p->n[1])*(p->n[0])*(p->ng[1]);
     else
       nmpibuffer=NDERV*(p->n[1])*(p->n[2])*(p->ng[1]);
   }
   else if((p->n[2])>=(p->n[0])  && (p->n[2])>=(p->n[1]))
   {
     if((p->n[0])>(p->n[1]))
       nmpibuffer=NDERV*(p->n[2])*(p->n[0])*(p->ng[2]);
     else
       nmpibuffer=NDERV*(p->n[2])*(p->n[1])*(p->ng[2]);
   }

#else
   if((p->n[0])>(p->n[1]))
    nmpibuffer=NDERV*(p->n[0])*(p->ng[0]);
   else
    nmpibuffer=NDERV*(p->n[1])*(p->ng[1]);
#endif
     gnmpirequest=0;
     gnmpibuffer=nmpibuffer;
     gnmpibuffermod=nmpibuffer*NVAR/NDERV;
     


     //gmpirequest=(MPI::Request *)calloc(numbuffers,sizeof(MPI::Request));

     gmpisendbuffer=(real *)calloc(nmpibuffer,sizeof(real));
     gmpirecvbuffer=(real *)calloc(nmpibuffer*numbuffers,sizeof(real));	
     

for(i=0;i<NDIM;i++)
{
     gmpisrcbufferl=(real **)calloc(NDIM,sizeof(real *));
     gmpisrcbufferr=(real **)calloc(NDIM,sizeof(real *));
     gmpitgtbufferl=(real **)calloc(NDIM,sizeof(real *));
     gmpitgtbufferr=(real **)calloc(NDIM,sizeof(real *));
}

for(i=0;i<NDIM;i++)
{
              switch(i)
              {
                 case 0:
#ifdef USE_SAC3D
     gmpisrcbufferl[i]=(real *)calloc( ((p->n[1])+2)*((p->n[2])+2)*(p->ng[0]),sizeof(real));
     gmpisrcbufferr[i]=(real *)calloc( ((p->n[1])+2)*((p->n[2])+2)*(p->ng[0]),sizeof(real ));
     gmpitgtbufferl[i]=(real *)calloc(((p->n[1])+2)*((p->n[2])+2)*(p->ng[0]),sizeof(real ));
     gmpitgtbufferr[i]=(real *)calloc(((p->n[1])+2)*((p->n[2])+2)*(p->ng[0]),sizeof(real ));
#else
     gmpisrcbufferl[i]=(real *)calloc(((p->n[1])+2)*(p->ng[0]),sizeof(real ));
     gmpisrcbufferr[i]=(real *)calloc(((p->n[1])+2)*(p->ng[0]),sizeof(real ));
     gmpitgtbufferl[i]=(real *)calloc(((p->n[1])+2)*(p->ng[0]),sizeof(real ));
     gmpitgtbufferr[i]=(real *)calloc(((p->n[1])+2)*(p->ng[0]),sizeof(real ));
#endif
                      
                      break;   
                 case 1:
#ifdef USE_SAC3D
     gmpisrcbufferl[i]=(real *)calloc(((p->n[0])+2)*((p->n[2])+2)*(p->ng[1]),sizeof(real ));
     gmpisrcbufferr[i]=(real *)calloc(((p->n[0])+2)*((p->n[2])+2)*(p->ng[1]),sizeof(real ));
     gmpitgtbufferl[i]=(real *)calloc(((p->n[0])+2)*((p->n[2])+2)*(p->ng[1]),sizeof(real ));
     gmpitgtbufferr[i]=(real *)calloc(((p->n[0])+2)*((p->n[2])+2)*(p->ng[1]),sizeof(real ));
#else
     gmpisrcbufferl[i]=(real *)calloc(((p->n[0])+2)*(p->ng[1]),sizeof(real ));
     gmpisrcbufferr[i]=(real *)calloc(((p->n[0])+2)*(p->ng[1]),sizeof(real ));
     gmpitgtbufferl[i]=(real *)calloc(((p->n[0])+2)*(p->ng[1]),sizeof(real ));
     gmpitgtbufferr[i]=(real *)calloc(((p->n[0])+2)*(p->ng[1]),sizeof(real ));
#endif
                      
                      break;
#ifdef USE_SAC3D         
                 case 2:
     gmpisrcbufferl[i]=(real *)calloc(((p->n[0])+2)*((p->n[1])+2)*(p->ng[2]),sizeof(real ));
     gmpisrcbufferr[i]=(real *)calloc(((p->n[0])+2)*((p->n[1])+2)*(p->ng[2]),sizeof(real ));
     gmpitgtbufferl[i]=(real *)calloc(((p->n[0])+2)*((p->n[1])+2)*(p->ng[2]),sizeof(real ));
     gmpitgtbufferr[i]=(real *)calloc(((p->n[0])+2)*((p->n[1])+2)*(p->ng[2]),sizeof(real ));
                      break;
#endif                             
                       }     
}    
     	

//comm.Barrier();


}



//!==============================================================================
//subroutine mpifinalize

//call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
//call MPI_FINALIZE(ierrmpi)

void mgpufinalize(struct params *p)
{
     //gwall_time = MPI_Wtime() - gwall_time;
     if ((p->ipe) == 0)
	  printf("\n Wall clock time = %f secs\n", gwall_time);
     //free(gmpisendbuffer);
     //free(gmpirecvbuffer);
     //free(gmpirequest);
     //MPI_Finalize();
}





//!==============================================================================
//subroutine mpisetnpeDipeD(name)

//! Set directional processor numbers and indexes based on a filename.
//! The filename contains _np followed by np^D written with 2 digit integers.
//! For example _np0203 means np1=2, np2=3 for 2D.

//! Extract and check the directional processor numbers and indexes
//! and concat the PE number to the input and output filenames

void mgpusetnpediped(struct params *p, char *string)
{    
  //we don't need to call this because the values p->pnpe[0],p->pnpe[1],p->pnpe[2]
  //have been set in the params file (they should be the same as the values in the filename)
  // for SAC these vaules are read from the filename
}


//!==============================================================================
//subroutine ipe2ipeD(qipe,qipe1,qipe2)
//
//! Convert serial processor index to directional processor indexes

//integer:: qipe1,qipe2, qipe
//!-----------------------------------------------------------------------------
//qipe1 = qipe - npe1*(qipe/npe1)
//qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 

void ipe2iped(struct params *p)
{

#ifdef USE_SAC_3D
//qipe1 = qipe - npe1*(qipe/npe1)
//qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 
//qipe3 = qipe/(npe1*npe2)
//(p->pipe[0])=(p->ipe)-(p->pnpe[0])*((p->ipe)/(p->pnpe[0]));
//(p->pipe[1])=((p->ipe)/(p->pnpe[0]))-(p->pnpe[1])*((p->ipe)/((p->pnpe[0])*(p->pnpe[1])));
//(p->pipe[2])=(p->ipe)/((p->pnpe[0])*(p->pnpe[1]));   

(p->pipe[2])=(p->ipe)/((p->pnpe[0])*(p->pnpe[1]));
(p->pipe[1])=((p->ipe)-((p->pipe[2])*(p->pnpe[0])*(p->pnpe[1])))/(p->pnpe[1]);
(p->pipe[0])=(p->ipe)-((p->pipe[2])*(p->pnpe[0])*(p->pnpe[1]))-((p->pipe[1])*(p->pnpe[1]));


//set upper boundary flags
//mpiupperB(1)=ipe1<npe1-1
//mpilowerB(1)=ipe1>0 

//mpiupperB(2)=ipe2<npe2-1
//mpilowerB(2)=ipe2>0 

(p->mpiupperb[0])=(p->pipe[0])<((p->pnpe[0])-1);
(p->mpiupperb[1])=(p->pipe[1])<((p->pnpe[1])-1);
(p->mpiupperb[2])=(p->pipe[2])<((p->pnpe[2])-1);

(p->mpilowerb[0])=(p->pipe[0])>0;
(p->mpilowerb[1])=(p->pipe[1])>0;
(p->mpilowerb[2])=(p->pipe[2])>0;
#else
//qipe1 = qipe - npe1*(qipe/npe1)
//qipe2 = qipe/npe1 - npe2*(qipe/(npe1*npe2)) 
//(p->pipe[0])=(p->ipe)-(p->pnpe[0])*(p->ipe)/(p->pnpe[0]);
//(p->pipe[1])=((p->ipe)/(p->pnpe[0]))-(p->pnpe[1])*(p->ipe)/((p->pnpe[0])*(p->pnpe[1]));
(p->pipe[1])=((p->ipe)/(p->pnpe[0]));
(p->pipe[0])=(p->ipe)-(p->pnpe[0])*(p->pipe[1]);

(p->mpiupperb[0])=(p->pipe[0])<((p->pnpe[0])-1);
(p->mpiupperb[1])=(p->pipe[1])<((p->pnpe[1])-1);

(p->mpilowerb[0])=(p->pipe[0])>0;
(p->mpilowerb[1])=(p->pipe[1])>0;


  if(p->ipe==0)
    for(int i=0; i<2;i++)
      printf("mpibc %d %d %d %d\n",i,p->pipe[i],p->mpiupperb[i],p->mpilowerb[i]);

#endif



    //ensure boundary set correctly 
    for(int ii=0; ii<NVAR; ii++)
    for(int idir=0; idir<NDIM; idir++)
    for(int ibound=0; ibound<2; ibound++)
    {
       if( ((p->boundtype[ii][idir][ibound])==0) && ((p->pnpe[idir])>0) )
                                     p->boundtype[ii][idir][ibound]=2;
       else if( (((p->mpiupperb[idir])==1) && (p->pipe[idir])<((p->pnpe[idir])-1) )  ||  ((p->mpiupperb[idir])!=1)  && ((p->pipe[idir])>0) )
                                     p->boundtype[ii][idir][ibound]=1;




      /* if( ((p->boundtype[ii][idir][ibound])==0) && ((p->pnpe[idir])>0) )
                                     p->boundtype[ii][idir][ibound]=2;
       else if( ((p->mpiupperb[idir])==1)    ||  ((p->mpilowerb[idir])==1) )
                                     p->boundtype[ii][idir][ibound]=1;*/



    }



}

//!==============================================================================
//subroutine ipeD2ipe(qipe1,qipe2,qipe)
//
//! Convert directional processor indexes to serial processor index

//include 'vacdef.f'

//integer:: qipe1,qipe2, qipe
//!-----------------------------------------------------------------------------
//qipe = qipe1  + npe1*qipe2

void iped2ipe(int *tpipe,int *tpnp, int *oipe)
{
  #ifdef USE_SAC_3D
  //qipe = qipe1  + npe1*qipe2  + npe1*npe2*qipe3
  (*oipe)=tpipe[0]+(tpnp[0])*(tpipe[1])+(tpnp[0])*(tpnp[1])*(tpipe[2]);
  #else
  (*oipe)=tpipe[0]+(tpnp[0])*(tpipe[1]);
  #endif

}

//!==============================================================================
//subroutine mpineighbors(idir,hpe,jpe)

//! Find the hpe and jpe processors on the left and right side of this processor 
//! in direction idir. The processor cube is taken to be periodic in every
//! direction.

//!-----------------------------------------------------------------------------
//hpe1=ipe1-kr(1,idir);hpe2=ipe2-kr(2,idir);
//jpe1=ipe1+kr(1,idir);jpe2=ipe2+kr(2,idir);

//if(hpe1<0)hpe1=npe1-1
//if(jpe1>=npe1)jpe1=0

//if(hpe2<0)hpe2=npe2-1
//if(jpe2>=npe2)jpe2=0

//call ipeD2ipe(hpe1,hpe2,hpe)
//call ipeD2ipe(jpe1,jpe2,jpe)

// Find the hpe and jpe processors on the left and right side of this processor 
// in direction idir. The processor cube is taken to be periodic in every
// direction.
void mgpuneighbours(int dir, struct params *p)
{
     int i;
     for(i=0; i<NDIM;i++)
     {
             
             (p->phpe[i])=(p->pipe[i])-(dir==i);
             (p->pjpe[i])=(p->pipe[i])+(dir==i);             
     }
     //printf("pcoords %d %d %d\n",p->ipe,p->pipe[0],p->pipe[1]);
     for(i=0; i<NDIM;i++)
     {
              if((p->phpe[i])<0) (p->phpe[i])=(p->pnpe[i])-1; 
              if((p->pjpe[i])<0) (p->pjpe[i])=(p->pnpe[i])-1; 
              if((p->phpe[i])>=(p->pnpe[i])) (p->phpe[i])=0; 
              if((p->pjpe[i])>=(p->pnpe[i])) (p->pjpe[i])=0;                     
     }
 // printf("lpcoords %d %d %d\n",p->ipe,p->phpe[0],p->phpe[1]);
//printf("rpcoords %d %d %d\n",p->ipe,p->pjpe[0],p->pjpe[1]);
   
     iped2ipe(p->phpe,p->pnpe,&(p->hpe));
     iped2ipe(p->pjpe,p->pnpe,&(p->jpe));
}

//!==============================================================================
//subroutine mpisend(nvar,var,ixmin1,ixmin2,ixmax1,ixmax2,qipe,iside)
//
//! Send var(ix^L,1:nvar) to processor qipe.
//! jside is 0 for min and 1 for max side of the grid for the sending PE
void mpisend(int nvar,real *var, int *ixmin, int *ixmax  ,int qipe,int iside, int dim, struct params *p)
{
    int n=0;
   int ivar,i1,i2,i3,bound;
   i3=0;

	switch(dim)
	{
		case 0:
		   for(ivar=0; ivar<nvar;ivar++)
		     for(i1=0;i1<=1;i1++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif

		      for(i2=0;i2<p->n[1];i2++)
		      {
			
                        bound=i1+2*(iside>0);
			 gmpisendbuffer[n]=var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)];



			//if((p->ipe==0) && ivar==rho && p->it != -1     /*&& iside==1 && (100*(p->ipe)+10*dim+iside)==101*/ )
			//{

                        //   bound=i1+2*(iside>0);
                        //    printf(" %d %d %d %lg  \n",bound,i2,i1,gmpisendbuffer[n]);

			//}
			//if((p->ipe==2) && ivar==pos2      /*&& iside==1 && (100*(p->ipe)+10*dim+iside)==101*/ )
                        //    printf(" %lg  \n",gmpisendbuffer[n]);


                         n++;






		      }
		   


		break;
		case 1:
		   for(ivar=0; ivar<nvar;ivar++)
                     for(i2=0;i2<=1;i2++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif

		     for(i1=0;i1<p->n[0];i1++)
		      
		      {
			
                        bound=i2+2*(iside>0);
			 gmpisendbuffer[n]=var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)];

			//if((p->ipe==0  ) && ivar==rho  && p->it != -1/* && ((p)->it)==2*/)
			//{
			 //for(int i=0;i<nvar;i++)
			 //  printf("mpiseend %d %d %d %lg \n",bound,i2,i1,gmpisendbuffer[n]);
                         //printf(" %d %d %d %lg ",bound,i2,i1,var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)]);
                         // ;//printf(" %d %d %d %lg  %lg\n",i1,i2,iside,var[sacencodempiw1 (p,i1, i2, i3, pos1,bound)],var[sacencodempiw1 (p,i1, i2, i3, pos2,bound)]);
			 //printf("\n");
			//}
                        n++;

		      }

 		/*for(i2=0;i2<=1;i2++)
                      for(i1=0;i1<p->n[0];i1++)
			if((p->ipe==2  )   && iside==1 && (100*(p->ipe)+10*dim+iside)==101 )
			{
			 //for(int i=0;i<nvar;i++)
                           bound=i2+2*(iside>0);
			   printf(" %d %d %d %lg  %lg\n",i1,i2,iside,var[sacencodempiw1 (p,i1, i2, i3, pos1,bound)],var[sacencodempiw1 (p,i1, i2, i3, pos2,bound)]);
			// printf("\n");
			}*/


		break;
		case 2:

		#ifdef USE_SAC3D
		   for(ivar=0; ivar<nvar;ivar++)
                     for(i3=0;i3<=1;i3++)
                     for(i2=0;p->n[1];i2++)
		     for(i1=0;i1<p->n[0];i1++)
		      		
						
		      {
			n++;
                        bound=i3+2*(iside>0);
			 gmpisendbuffer[n]=var[sacencodempiw2 (p,i1, i2, i3, ivar,bound)];
		      }




		#endif

		      
		break;
	}






//if(p->ipe==1  && dim==0)
//      printf("ipe %d send tag %d nb %d  to %d %d %d\n",p->ipe,100*((p->ipe)+1)+10*(dim+1)+(iside==0?1:0),n,qipe,iside,dim);
   
 //  comm.Rsend(gmpisendbuffer, n, MPI_DOUBLE_PRECISION, qipe, 100*((p->ipe)+1)+10*(dim+1)+(iside==0?1:0));


}





void mpisendmod(int nvar,real *var, int *ixmin, int *ixmax  ,int qipe,int iside, int dim, struct params *p)
{
    int n=0;
   int ivar,i1,i2,i3,bound;
   i3=0;

	switch(dim)
	{
		case 0:
		   for(ivar=0; ivar<nvar;ivar++)
		     for(i1=0;i1<=1;i1++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif

		      for(i2=0;i2<p->n[1];i2++)
		      {
			
                        bound=i1+2*(iside>0);
			 gmpisendbuffer[n]=var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)];



			//if((p->ipe==0) && ivar==rho /*&& p->it != -1  */   && iside==1/* && (100*(p->ipe)+10*dim+iside)==101*/ )
			//{

                        //   bound=i1+2*(iside>0);
                        //    printf(" %d %d %d %d %d %lg  \n",bound,i2,i1,n,sacencodempiw0 (p,i1, i2, i3, ivar,bound),gmpisendbuffer[n]);

                        // printf("3 %d %d %d %lg  %lg\n",i1,i2,n,var[sacencodempiw0 (p,i1, i2, i3, rho,bound)],gmpisendbuffer[n]);


			//}
			//if((p->ipe==2) && ivar==pos2      /*&& iside==1 && (100*(p->ipe)+10*dim+iside)==101*/ )
                        //    printf(" %lg  \n",gmpisendbuffer[n]);


                         n++;






		      }
		   


		break;
		case 1:
		   for(ivar=0; ivar<nvar;ivar++)
                     for(i2=0;i2<=1;i2++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif

		     for(i1=0;i1<p->n[0];i1++)
		      
		      {
			
                        bound=i2+2*(iside>0);
			 gmpisendbuffer[n]=var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)];

			//if((p->ipe==0  ) && ivar==rho  /*&& p->it != -1/* && ((p)->it)==2*/)
			//{
			 //for(int i=0;i<nvar;i++)
			 //  printf("mpiseend %d %d %d %lg \n",bound,i2,i1,gmpisendbuffer[n]);
                         //printf(" %d %d %d %lg ",bound,i2,i1,var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)]);
                         //printf(" %d %d %d %lg  %lg\n",i1,i2,iside,var[sacencodempiw1 (p,i1, i2, i3, rho,bound)],var[sacencodempiw1 (p,i1, i2, i3, rho,bound)]);
			 //printf("\n");
			//}
                        n++;

		      }

 		/*for(i2=0;i2<=1;i2++)
                      for(i1=0;i1<p->n[0];i1++)
			if((p->ipe==2  )   && iside==1 && (100*(p->ipe)+10*dim+iside)==101 )
			{
			 //for(int i=0;i<nvar;i++)
                           bound=i2+2*(iside>0);
			   printf(" %d %d %d %lg  %lg\n",i1,i2,iside,var[sacencodempiw1 (p,i1, i2, i3, pos1,bound)],var[sacencodempiw1 (p,i1, i2, i3, pos2,bound)]);
			// printf("\n");
			}*/


		break;
		case 2:

		#ifdef USE_SAC3D
		   for(ivar=0; ivar<nvar;ivar++)
                     for(i3=0;i3<=1;i3++)
                     for(i2=0;p->n[1];i2++)
		     for(i1=0;i1<p->n[0];i1++)
		      		
						
		      {
			n++;
                        bound=i3+2*(iside>0);
			 gmpisendbuffer[n]=var[sacencodempiw2 (p,i1, i2, i3, ivar,bound)];
		      }


		#endif

		      
		break;
	}






/*if(p->ipe==0  && dim==0 && iside==0)
{
  for(i1=0; i1<5120;i1++)
      printf("ipe %d send tag %d nb %d  to %d %d %d %lg %lg\n",p->ipe,100*((p->ipe)+1)+10*(dim+1)+(iside==0?1:0),n,qipe,iside,i1,gmpisendbuffer[i1],var[i1]);
  printf("end\n\n");
}*/
   
 //  comm.Rsend(gmpisendbuffer, n, MPI_DOUBLE_PRECISION, qipe, 100*((p->ipe)+1)+10*(dim+1)+(iside==0?1:0));

//if((p->ipe==1  )   )
//			printf("%d %d %d\n",p->ipe,iside,n);
}



//!==============================================================================
//subroutine mpirecvbuffer(nvar,ixmin1,ixmin2,ixmax1,ixmax2,qipe,iside)
//
//! receive buffer for a ghost cell region of size ix^L sent from processor qipe
//! and sent from side iside of the grid
//
//include 'vacdef.f'
//
//integer:: nvar, ixmin1,ixmin2,ixmax1,ixmax2, qipe, iside, n
//
//integer :: nmpirequest, mpirequests(2)
//integer :: mpistatus(MPI_STATUS_SIZE,2)
//common /mpirecv/ nmpirequest,mpirequests,mpistatus
//!----------------------------------------------------------------------------
void mpirecvbuffer(int nvar,int *ixmin, int *ixmax  ,int qipe,int iside,int dim, struct params *p)
{
int nrecv;
/*#ifdef USE_SAC3D
   nrecv = nvar* (ixmax[0]-ixmin[0]+1)*(ixmax[1]-ixmin[1]+1)*(ixmax[2]-ixmin[2]+1);
#else
   nrecv = nvar* (ixmax[0]-ixmin[0]+1)*(ixmax[1]-ixmin[1]+1);
#endif*/

	switch(dim)
	{
		case 0:
			#ifdef USE_SAC3D
			   nrecv = 2*nvar* (p->n[1])*(p->n[2]);
			#else
			   nrecv = 2*nvar* (p->n[1]);
			#endif
		break;
		case 1:
			#ifdef USE_SAC3D
			   nrecv = 2*nvar* (p->n[0])*(p->n[2]);
			#else
			   nrecv = 2*nvar* (p->n[0]);
			#endif
		break;
		case 2:
			#ifdef USE_SAC3D
			   nrecv = 2*nvar* (p->n[1])*(p->n[0]);
			#endif
		break;
	}


gnmpirequest++;
//  if((p->ipe)==0  && dim==0)
//      printf("ipe %d recv tag %d nb %d  to %d  %d %d\n",p->ipe, 100*(qipe+1)+10*(dim+1)+iside/*(iside==0?1:0)*/ ,nrecv,qipe,iside,dim);

//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer+(iside*gnmpibuffer),nrecv,MPI_DOUBLE_PRECISION,qipe,10*(p->ipe)+iside);
//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer+(iside*gnmpibuffer),nrecv,MPI_DOUBLE_PRECISION,qipe,MPI_ANY_TAG);
//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer+(2*iside*gnmpibuffer),nrecv,MPI_DOUBLE_PRECISION,qipe,100*(qipe+1)+10*(dim+1)+iside/**(iside==0?1:0)*/);
//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer,nrecv,MPI_DOUBLE_PRECISION,qipe,100*(qipe+1)+10*(dim+1)+iside/**(iside==0?1:0)*/);
}





void mpirecvbuffermod(int nvar,int *ixmin, int *ixmax  ,int qipe,int iside,int dim, struct params *p)
{
int nrecv;
/*#ifdef USE_SAC3D
   nrecv = nvar* (ixmax[0]-ixmin[0]+1)*(ixmax[1]-ixmin[1]+1)*(ixmax[2]-ixmin[2]+1);
#else
   nrecv = nvar* (ixmax[0]-ixmin[0]+1)*(ixmax[1]-ixmin[1]+1);
#endif*/

	switch(dim)
	{
		case 0:
			#ifdef USE_SAC3D
			   nrecv = 2*nvar* (p->n[1])*(p->n[2]);
			#else
			   nrecv = 2*nvar* (p->n[1]);
			#endif
		break;
		case 1:
			#ifdef USE_SAC3D
			   nrecv = 2*nvar* (p->n[0])*(p->n[2]);
			#else
			   nrecv = 2*nvar* (p->n[0]);
			#endif
		break;
		case 2:
			#ifdef USE_SAC3D
			   nrecv = 2*nvar* (p->n[1])*(p->n[0]);
			#endif
		break;
	}


gnmpirequest++;
//  if((p->ipe)==0  && dim==0)
//      printf("ipe %d recv tag %d nb %d  to %d  %d %d\n",p->ipe, 100*(qipe+1)+10*(dim+1)+iside/*(iside==0?1:0)*/ ,nrecv,qipe,iside,dim);

//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer+(iside*gnmpibuffer),nrecv,MPI_DOUBLE_PRECISION,qipe,10*(p->ipe)+iside);
//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer+(iside*gnmpibuffer),nrecv,MPI_DOUBLE_PRECISION,qipe,MPI_ANY_TAG);
//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer+(2*iside*gnmpibuffermod),nrecv,MPI_DOUBLE_PRECISION,qipe,100*(qipe+1)+10*(dim+1)+iside/**(iside==0?1:0)*/);
//gmpirequest[gnmpirequest]=comm.Irecv(gmpirecvbuffer,nrecv,MPI_DOUBLE_PRECISION,qipe,100*(qipe+1)+10*(dim+1)+iside/**(iside==0?1:0)*/);
}




//!==============================================================================
//subroutine mpibuffer2var(iside,nvar,var,ixmin1,ixmin2,ixmax1,ixmax2)
//
//! Copy mpibuffer(:,iside) into var(ix^L,1:nvar)
//include 'vacdef.f'
//
//integer :: nvar
//double precision:: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nvar)
//integer:: ixmin1,ixmin2,ixmax1,ixmax2,iside,n,ix1,ix2,ivar
//!-----------------------------------------------------------------------------
void mpibuffer2var(int iside,int nvar,real *var, int *ixmin, int *ixmax, int dim, struct params *p)
{
   int n=0;
   int ivar,i1,i2,i3,bound;
   i3=0;

	switch(dim)
	{
		case 0:
		   for(ivar=0; ivar<nvar;ivar++)
		     for(i1=0;i1<=1;i1++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif
		      for(i2=0;i2<p->n[1];i2++)

		      {
			
                        //bound=i1+iside+2*(iside>0);
                        bound=i1+2*(iside==0?1:0);
			// var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n/*+iside*gnmpibuffer*/];
 			var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer];
                        //n++;
                       //iside=1;
                       //if(/*iside==0  && p->it != -1 &&*/ p->ipe==0  && ivar==rho /*&& (100+10*dim+(iside==0?1:0))==101*/)
                       //   printf("mpib2var %d %d %d %lg\n",i1 ,i2 , iside,gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer]);
                       n++;
                         
			//printf("\n");

		      }
                      n=0;
                   /* if(p->ipe==2)
			// for(i1=0;i1<=1;i1++)
			 //    for(i2=0;i2<p->n[1];i2++)
                            for(i1=0; i1<4*gnmpibuffer; i1++)
                             {
                               
                                 printf("%d %d %d %lg\n",i1 ,i2 , iside,gmpirecvbuffer[i1]);
                                 n++;
                             }*/
		break;
		case 1:
		   for(ivar=0; ivar<nvar;ivar++)
                    for(i2=0;i2<=1;i2++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif
		     for(i1=0;i1<p->n[0];i1++)
		      

		      {
			
                        bound=i2+2*(iside==0?1:0);
			 var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer];


                      // if(/*iside==0  && p->it != -1 &&*/ p->ipe==0  && ivar==rho /*&& (100+10*dim+(iside==0?1:0))==101*/)
                      //    printf("mpib2var %d %d %d %lg\n",i1 ,i2 , iside,var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)]);

                       // if(/*iside==0  &&*/ p->ipe==2  && ivar==pos1 || ivar==pos2 /*&& (100+10*dim+(iside==0?1:0))==101*/)
                       //   printf("v %d %d %d %d %lg\n",ivar,i1 ,i2 , iside,gmpirecvbuffer[n+iside*gnmpibuffer]);

                        n++;

		      }

		break;
		case 2:

		#ifdef USE_SAC3D
		   for(ivar=0; ivar<nvar;ivar++)
			for(i3=0;i3<=1;i3++)			
		      for(i2=0;p->n[1];i2++)		
		     for(i1=0;i1<p->n[0];i1++)
		      {
			
                        bound=i3+2*(iside==0?1:0);
			 var[sacencodempiw2 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer];
                          n++;
		      }


		#endif

		      
		break;
	}
}



void mpibuffer2varmod(int iside,int nvar,real *var, int *ixmin, int *ixmax, int dim, struct params *p)
{
   int n=0;
   int ivar,i1,i2,i3,bound;
   i3=0;
                 //  ivar=0;
                //  for(n=0;n<gnmpibuffermod;n++)
                 //     if(dim==0  &&/* p->it != -1 &&*/ p->ipe==0  /*&& ivar==rho && (100+10*dim+(iside==0?1:0))==101*/){
                  //         printf("mpib2var %d %d %lg\n" , iside,n,gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffermod]);
			// printf("mpib2var %d %d %d %lg\n",i1 ,i2 , iside,var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]);
                      //   printf("\n");
                   //    }
    

         n=0;

	switch(dim)
	{
		case 0:
		   for(ivar=0; ivar<nvar;ivar++)
		     for(i1=0;i1<=1;i1++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif
		      for(i2=0;i2<p->n[1];i2++)

		      {
			
                        //bound=i1+iside+2*(iside>0);
                        bound=i1+2*(iside==0?1:0);
			// var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n/*+iside*gnmpibuffer*/];
 			//var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+(iside==0?1:0)*gnmpibuffer/nvar];
                       // var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer];
                        var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffermod];
                        //n++;
                       //iside=1;
                      // if(/*iside==0  && p->it != -1 &&*/ p->ipe==0  && ivar==rho /*&& (100+10*dim+(iside==0?1:0))==101*/){
                      //    printf("mpib2vart %d %d %d %d %lg\n",i1 ,i2 , iside,n,gmpirecvbuffer[n+2*iside*gnmpibuffer]);
			// printf("mpib2var %d %d %d %lg\n",i1 ,i2 , iside,var[sacencodempiw0 (p,i1, i2, i3, ivar,bound)]);
                      //   printf("\n");
                     //  }
                       n++;
                         
			

		      }
                     // n=0;
                   /* if(p->ipe==0)
			// for(i1=0;i1<=1;i1++)
			 //    for(i2=0;i2<p->n[1];i2++)
                            for(i1=0; i1<4*gnmpibuffer; i1++)
                             {

                               
                                 printf("%d %d %d %lg\n",i1 ,i2 , iside,gmpirecvbuffer[i1]);
                                 n++;
                             }*/
		break;
		case 1:
		   for(ivar=0; ivar<nvar;ivar++)
                    for(i2=0;i2<=1;i2++)
		#ifdef USE_SAC3D
			for(i3=0;i3<p->n[2];i3++)
		#endif
		     for(i1=0;i1<p->n[0];i1++)
		      

		      {
			
                        bound=i2+2*(iside==0?1:0);
			// var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer];
			var[sacencodempiw1 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffermod];


                       //if(/*iside==0  && p->it != -1 &&*/ p->ipe==0  && ivar==rho /*&& (100+10*dim+(iside==0?1:0))==101*/){
                       //   printf("mpib2var %d %d %d %lg\n",i1 ,i2 , iside,gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer]);
                       //   printf("\n");
                       // }




                       // if(/*iside==0  &&*/ p->ipe==2  && ivar==pos1 || ivar==pos2 /*&& (100+10*dim+(iside==0?1:0))==101*/)
                       //   printf("v %d %d %d %d %lg\n",ivar,i1 ,i2 , iside,gmpirecvbuffer[n+iside*gnmpibuffer]);

                        n++;

		      }

		break;
		case 2:

		#ifdef USE_SAC3D
		   for(ivar=0; ivar<nvar;ivar++)
			for(i3=0;i3<=1;i3++)			
		      for(i2=0;p->n[1];i2++)		
		     for(i1=0;i1<p->n[0];i1++)
		      {
			
                        bound=i3+2*(iside==0?1:0);
			// var[sacencodempiw2 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffer];
 			var[sacencodempiw2 (p,i1, i2, i3, ivar,bound)]=gmpirecvbuffer[n+2*(iside==0?1:0)*gnmpibuffermod];
                          n++;
		      }


		#endif

		      
		break;
	}
}




void gpusync()
{
//printf("mpisync\n");
//comm.Barrier();
;
}









//!==============================================================================
//subroutine mpibound(nvar,var)
//
//! Fill in ghost cells of var(ixG,nvar) from other processors
//
//include 'vacdef.f'
//
//integer :: nvar
//double precision :: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nvar)
//
//! processor indexes for left and right neighbors
//integer :: hpe,jpe
//! index limits for the left and right side mesh and ghost cells 
//integer :: ixLMmin1,ixLMmin2,ixLMmax1,ixLMmax2, ixRMmin1,ixRMmin2,ixRMmax1,&
//   ixRMmax2, ixLGmin1,ixLGmin2,ixLGmax1,ixLGmax2, ixRGmin1,ixRGmin2,ixRGmax1,&
//   ixRGmax2
//logical :: periodic


void mpibound(int nvar,  real *var1, real *var2, real *var3, struct params *p, int idir)
{
   int i;

   int ixlgmin[NDIM],ixlgmax[NDIM];
   int ixrgmin[NDIM],ixrgmax[NDIM];

   int ixlmmin[NDIM],ixlmmax[NDIM];
   int ixrmmin[NDIM],ixrmmax[NDIM];


if((p->pnpe[0])>1   && idir==0)
{
   gnmpirequest=0;
   for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   //periodic=typeB(1,2*1)=='mpiperiod'
   //! Left and right side ghost cell regions (target)
   for(i=0;i<NDIM;i++)
   {
	ixlgmin[i]=0;
	ixlgmax[i]=(p->n[i])-1;
	ixrgmin[i]=0;
	ixrgmax[i]=(p->n[i])-1;
   }
      ixlgmax[0]=1;
   ixrgmin[0]=(p->n[0])-2;
   //! Left and right side mesh cell regions (source)
   for(i=0;i<NDIM;i++)
   {
	ixlmmin[i]=0;
	ixlmmax[i]=(p->n[i])-1;
	ixrmmin[i]=0;
	ixrmmax[i]=(p->n[i])-1;
   }
   ixlmmin[0]=2;   
   ixlmmax[0]=1;
   ixrmmax[0]=(p->n[0])-3;
   ixrmmin[0]=(p->n[0])-4; 
     //! Obtain left and right neighbor processors for this direction
   //call mpineighbors(1,hpe,jpe)
   
   mgpuneighbours(0, p);
   //already computed initially use phpe pjpe arrays

   //! receive right (2) boundary from left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,hpe,2)*/

   //if(p->ipe==0)
   //  printf("ipe %d  recv right (from left) %d recv left (from right neigh) %d\n",p->ipe,p->hpe,p->jpe);
   if(((p->mpilowerb[0])==1) ||  /*((p->boundtype[0][0][0])==0)||*/  ((p->boundtype[0][0][0])==2))
	mpirecvbuffer(nvar,ixrmmin,ixrmmax,p->hpe,0,0,p);

   //! receive left (1) boundary from right neighbor jpe
   //if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,jpe,1)
   if(((p->mpiupperb[0])==1) || /* ((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
             mpirecvbuffer(nvar,ixlmmin,ixlmmax,p->jpe,1,0,p);


  
   //! Wait for all receives to be posted
   //call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
//printf("ipe %d barrier before send\n",p->ipe);
  // comm.Barrier();
   
   //! Ready send left (1) boundary to left neighbor hpe

   //if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,hpe,1)
   if(((p->mpilowerb[0])==1) ||  /*((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
            mpisend(nvar,var1,ixlmmin,ixlmmax,p->hpe,0,0,p);

   //! Ready send right (2) boundary to right neighbor
   //if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,jpe,2)
   if(((p->mpiupperb[0])==1) || /* ((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
             mpisend(nvar,var1,ixrmmin,ixrmmax,p->jpe,1,0,p);

   //! Wait for messages to arrive
   //call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
  
   /*request.Waitall(gnmpirequest,gmpirequest);*/

   //! Copy buffer received from right (2) physical cells into left ghost cells
   //if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
   //   ixLGmin2,ixLGmax1,ixLGmax2)
   if(((p->mpilowerb[0])==1) ||  /*((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
             mpibuffer2var(1,nvar,var1,ixlgmin,ixlgmax,0,p);
   //! Copy buffer received from left (1) physical cells into right ghost cells
   //if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
   //   ixRGmin2,ixRGmax1,ixRGmax2)
   if(((p->mpiupperb[0])==1) || /* ((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
             mpibuffer2var(0,nvar,var1,ixrgmin,ixrgmax,0,p);    
  
 

}



 /*  comm.Barrier();  */


//printf("to here1 %d\n");
if((p->pnpe[1])>1   && idir==1)
{
  gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   //periodic=typeB(1,2*1)=='mpiperiod'
   //! Left and right side ghost cell regions (target)
   for(i=0;i<NDIM;i++)
   {
	ixlgmin[i]=0;
	ixlgmax[i]=(p->n[i])-1;
	ixrgmin[i]=0;
	ixrgmax[i]=(p->n[i])-1;
   }
   ixlgmax[1]=1;
   ixrgmin[1]=(p->n[1])-2;
   //! Left and right side mesh cell regions (source)
   for(i=0;i<NDIM;i++)
   {
	ixlmmin[i]=0;
	ixlmmax[i]=(p->n[i])-1;
	ixrmmin[i]=0;
	ixrmmax[i]=(p->n[i])-1;
   }
   ixlmmin[1]=2;   
   ixlmmax[1]=1;
   ixrmmax[1]=(p->n[1])-3;
   ixrmmin[1]=(p->n[1])-4; 
     //! Obtain left and right neighbor processors for this direction
   //call mpineighbors(1,hpe,jpe)
   mgpuneighbours(1, p);

   //already computed initially use phpe pjpe arrays

   //! receive right (2) boundary from left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,hpe,2)
   if(((p->mpilowerb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpirecvbuffer(nvar,ixrmmin,ixrmmax,p->hpe,0,1,p);
   //! receive left (1) boundary from right neighbor jpe
   //if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,jpe,1)
   if(((p->mpiupperb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpirecvbuffer(nvar,ixlmmin,ixlmmax,p->jpe,1,1,p);
   //! Wait for all receives to be posted
   //call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   /*comm.Barrier();*/
   //! Ready send left (1) boundary to left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,hpe,1)
   if(((p->mpilowerb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpisend(nvar,var2,ixlmmin,ixlmmax,p->hpe,0,1,p);
   //! Ready send right (2) boundary to right neighbor
   //if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,jpe,2)
   if(((p->mpiupperb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpisend(nvar,var2,ixrmmin,ixrmmax,p->jpe,1,1,p);
   //! Wait for messages to arrive
   //call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   /*request.Waitall(gnmpirequest,gmpirequest);*/

   //! Copy buffer received from right (2) physical cells into left ghost cells
   //if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
   //   ixLGmin2,ixLGmax1,ixLGmax2)
   if(((p->mpilowerb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpibuffer2var(1,nvar,var2,ixlgmin,ixlgmax,1,p);
   //! Copy buffer received from left (1) physical cells into right ghost cells
   //if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
   //   ixRGmin2,ixRGmax1,ixRGmax2)
   if(((p->mpiupperb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpibuffer2var(0,nvar,var2,ixrgmin,ixrgmax,1,p);
}

 /*comm.Barrier();*/




#ifdef USE_SAC3D
 /*comm.Barrier();*/
if((p->pnpe[2])>1)
{
  gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   //periodic=typeB(1,2*1)=='mpiperiod'
   //! Left and right side ghost cell regions (target)
   for(i=0;i<NDIM;i++)
   {
	ixlgmin[i]=0;
	ixlgmax[i]=(p->n[i])-1;
	ixrgmin[i]=0;
	ixrgmax[i]=(p->n[i])-1;
   }
   ixlgmax[2]=1;
   ixrgmin[2]=(p->n[2])-2;
   //! Left and right side mesh cell regions (source)
   for(i=0;i<NDIM;i++)
   {
	ixlmmin[i]=0;
	ixlmmax[i]=(p->n[i])-1;
	ixrmmin[i]=0;
	ixrmmax[i]=(p->n[i])-1;
   }
   ixlmmin[2]=2;   
   ixlmmax[2]=1;
   ixrmmax[2]=(p->n[2])-3;
   ixrmmin[2]=(p->n[2])-4; 
     //! Obtain left and right neighbor processors for this direction
   //call mpineighbors(1,hpe,jpe)
   mgpuneighbours(2, p);

   //already computed initially use phpe pjpe arrays

   //! receive right (2) boundary from left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,hpe,2)
   if(((p->mpilowerb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpirecvbuffer(nvar,ixrmmin,ixrmmax,p->hpe,0,2,p);
   //! receive left (1) boundary from right neighbor jpe
   //if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,jpe,1)
   if(((p->mpiupperb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpirecvbuffer(nvar,ixlmmin,ixlmmax,p->jpe,1,2,p);
   //! Wait for all receives to be posted
   //call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   /*comm.Barrier();*/
   //! Ready send left (1) boundary to left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,hpe,1)
   if(((p->mpilowerb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpisend(nvar,var3,ixlmmin,ixlmmax,p->hpe,0,2,p);
   //! Ready send right (2) boundary to right neighbor
   //if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,jpe,2)
   if(((p->mpiupperb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpisend(nvar,var3,ixrmmin,ixrmmax,p->jpe,1,2,p);
   //! Wait for messages to arrive
   //call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   /*request.Waitall(gnmpirequest,gmpirequest);*/

   //! Copy buffer received from right (2) physical cells into left ghost cells
   //if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
   //   ixLGmin2,ixLGmax1,ixLGmax2)
   if(((p->mpilowerb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpibuffer2var(1,nvar,var3,ixlgmin,ixlgmax,2,p);
   //! Copy buffer received from left (1) physical cells into right ghost cells
   //if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
   //   ixRGmin2,ixRGmax1,ixRGmax2)
   if(((p->mpiupperb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpibuffer2var(0,nvar,var3,ixrgmin,ixrgmax,2,p);
}


#endif


}






void mpiboundmod(int nvar,  real *var1, real *var2, real *var3, struct params *p, int idir)
{
   int i;

   int ixlgmin[NDIM],ixlgmax[NDIM];
   int ixrgmin[NDIM],ixrgmax[NDIM];

   int ixlmmin[NDIM],ixlmmax[NDIM];
   int ixrmmin[NDIM],ixrmmax[NDIM];


if((p->pnpe[0])>1   && idir==0)
{
   gnmpirequest=0;
   for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   //periodic=typeB(1,2*1)=='mpiperiod'
   //! Left and right side ghost cell regions (target)
   for(i=0;i<NDIM;i++)
   {
	ixlgmin[i]=0;
	ixlgmax[i]=(p->n[i])-1;
	ixrgmin[i]=0;
	ixrgmax[i]=(p->n[i])-1;
   }
      ixlgmax[0]=1;
   ixrgmin[0]=(p->n[0])-2;
   //! Left and right side mesh cell regions (source)
   for(i=0;i<NDIM;i++)
   {
	ixlmmin[i]=0;
	ixlmmax[i]=(p->n[i])-1;
	ixrmmin[i]=0;
	ixrmmax[i]=(p->n[i])-1;
   }
   ixlmmin[0]=2;   
   ixlmmax[0]=1;
   ixrmmax[0]=(p->n[0])-3;
   ixrmmin[0]=(p->n[0])-4; 
     //! Obtain left and right neighbor processors for this direction
   //call mpineighbors(1,hpe,jpe)
   
   mgpuneighbours(0, p);
   //already computed initially use phpe pjpe arrays

   //! receive right (2) boundary from left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,hpe,2)*/

   //if(p->ipe==0)
   //  printf("ipe %d  recv right (from left) %d recv left (from right neigh) %d\n",p->ipe,p->hpe,p->jpe);
   if(((p->mpilowerb[0])==1) ||  /*((p->boundtype[0][0][0])==0)||*/  ((p->boundtype[0][0][0])==2))
	mpirecvbuffermod(nvar,ixrmmin,ixrmmax,p->hpe,0,0,p);

   //! receive left (1) boundary from right neighbor jpe
   //if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,jpe,1)
   if(((p->mpiupperb[0])==1) || /* ((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
             mpirecvbuffermod(nvar,ixlmmin,ixlmmax,p->jpe,1,0,p);


  
   //! Wait for all receives to be posted
   //call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
//printf("ipe %d barrier before send\n",p->ipe);
  // comm.Barrier();
   
   //! Ready send left (1) boundary to left neighbor hpe

   //if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,hpe,1)
   if(((p->mpilowerb[0])==1) ||  /*((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
            mpisendmod(nvar,var1,ixlmmin,ixlmmax,p->hpe,0,0,p);

   //! Ready send right (2) boundary to right neighbor
   //if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,jpe,2)
   if(((p->mpiupperb[0])==1) || /* ((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
             mpisendmod(nvar,var1,ixrmmin,ixrmmax,p->jpe,1,0,p);

   //! Wait for messages to arrive
   //call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
  
  // request.Waitall(gnmpirequest,gmpirequest);

   //! Copy buffer received from right (2) physical cells into left ghost cells
   //if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
   //   ixLGmin2,ixLGmax1,ixLGmax2)
   
//comm.Barrier();
if(((p->mpilowerb[0])==1) ||  /*((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
   {
     //if((p->ipe)==0)
      //  printf("lowerb");
             mpibuffer2varmod(1,nvar,var1,ixlgmin,ixlgmax,0,p);
   }
   //! Copy buffer received from left (1) physical cells into right ghost cells
   //if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
   //   ixRGmin2,ixRGmax1,ixRGmax2)
//comm.Barrier();
   if(((p->mpiupperb[0])==1) || /* ((p->boundtype[0][0][0])==0)|| */ ((p->boundtype[0][0][0])==2))
   {
     //if((p->ipe)==0)
     //   printf("upperb");
             mpibuffer2varmod(0,nvar,var1,ixrgmin,ixrgmax,0,p);  
   }  
  
 

}



   /*comm.Barrier();*/


//printf("to here1 %d\n");
if((p->pnpe[1])>1   && idir==1)
{
  gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   //periodic=typeB(1,2*1)=='mpiperiod'
   //! Left and right side ghost cell regions (target)
   for(i=0;i<NDIM;i++)
   {
	ixlgmin[i]=0;
	ixlgmax[i]=(p->n[i])-1;
	ixrgmin[i]=0;
	ixrgmax[i]=(p->n[i])-1;
   }
   ixlgmax[1]=1;
   ixrgmin[1]=(p->n[1])-2;
   //! Left and right side mesh cell regions (source)
   for(i=0;i<NDIM;i++)
   {
	ixlmmin[i]=0;
	ixlmmax[i]=(p->n[i])-1;
	ixrmmin[i]=0;
	ixrmmax[i]=(p->n[i])-1;
   }
   ixlmmin[1]=2;   
   ixlmmax[1]=1;
   ixrmmax[1]=(p->n[1])-3;
   ixrmmin[1]=(p->n[1])-4; 
     //! Obtain left and right neighbor processors for this direction
   //call mpineighbors(1,hpe,jpe)
   mgpuneighbours(1, p);

   //already computed initially use phpe pjpe arrays

   //! receive right (2) boundary from left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,hpe,2)
   if(((p->mpilowerb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpirecvbuffermod(nvar,ixrmmin,ixrmmax,p->hpe,0,1,p);
   //! receive left (1) boundary from right neighbor jpe
   //if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,jpe,1)
   if(((p->mpiupperb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpirecvbuffermod(nvar,ixlmmin,ixlmmax,p->jpe,1,1,p);
   //! Wait for all receives to be posted
   //call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   /*comm.Barrier();*/
   //! Ready send left (1) boundary to left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,hpe,1)
   if(((p->mpilowerb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpisendmod(nvar,var2,ixlmmin,ixlmmax,p->hpe,0,1,p);
   //! Ready send right (2) boundary to right neighbor
   //if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,jpe,2)
   if(((p->mpiupperb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpisendmod(nvar,var2,ixrmmin,ixrmmax,p->jpe,1,1,p);
   //! Wait for messages to arrive
   //call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   /*request.Waitall(gnmpirequest,gmpirequest);*/

   //! Copy buffer received from right (2) physical cells into left ghost cells
   //if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
   //   ixLGmin2,ixLGmax1,ixLGmax2)
   if(((p->mpilowerb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpibuffer2varmod(1,nvar,var2,ixlgmin,ixlgmax,1,p);
   //! Copy buffer received from left (1) physical cells into right ghost cells
   //if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
   //   ixRGmin2,ixRGmax1,ixRGmax2)
   if(((p->mpiupperb[1])==1) ||  /*((p->boundtype[0][1][0])==0)||*/  ((p->boundtype[0][0][0])==2))
             mpibuffer2varmod(0,nvar,var2,ixrgmin,ixrgmax,1,p);
}

 /*comm.Barrier();*/




#ifdef USE_SAC3D
 /*comm.Barrier();*/
if((p->pnpe[2])>1)
{
  gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   //periodic=typeB(1,2*1)=='mpiperiod'
   //! Left and right side ghost cell regions (target)
   for(i=0;i<NDIM;i++)
   {
	ixlgmin[i]=0;

	ixlgmax[i]=(p->n[i])-1;
	ixrgmin[i]=0;
	ixrgmax[i]=(p->n[i])-1;
   }
   ixlgmax[2]=1;
   ixrgmin[2]=(p->n[2])-2;
   //! Left and right side mesh cell regions (source)
   for(i=0;i<NDIM;i++)
   {
	ixlmmin[i]=0;
	ixlmmax[i]=(p->n[i])-1;
	ixrmmin[i]=0;
	ixrmmax[i]=(p->n[i])-1;
   }
   ixlmmin[2]=2;   
   ixlmmax[2]=1;
   ixrmmax[2]=(p->n[2])-3;
   ixrmmin[2]=(p->n[2])-4; 
     //! Obtain left and right neighbor processors for this direction
   //call mpineighbors(1,hpe,jpe)
   mgpuneighbours(2, p);

   //already computed initially use phpe pjpe arrays

   //! receive right (2) boundary from left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpirecvbuffer(nvar,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,hpe,2)
   if(((p->mpilowerb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpirecvbuffermod(nvar,ixrmmin,ixrmmax,p->hpe,0,2,p);
   //! receive left (1) boundary from right neighbor jpe
   //if(mpiupperB(1) .or. periodic)call mpirecvbuffer(nvar,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,jpe,1)
   if(((p->mpiupperb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpirecvbuffermod(nvar,ixlmmin,ixlmmax,p->jpe,1,2,p);
   //! Wait for all receives to be posted
   //call MPI_BARRIER(MPI_COMM_WORLD,ierrmpi)
   /*comm.Barrier();*/
   //! Ready send left (1) boundary to left neighbor hpe
   //if(mpilowerB(1) .or. periodic)call mpisend(nvar,var,ixLMmin1,ixLMmin2,&
   //   ixLMmax1,ixLMmax2,hpe,1)
   if(((p->mpilowerb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpisendmod(nvar,var3,ixlmmin,ixlmmax,p->hpe,0,2,p);
   //! Ready send right (2) boundary to right neighbor
   //if(mpiupperB(1) .or. periodic)call mpisend(nvar,var,ixRMmin1,ixRMmin2,&
   //   ixRMmax1,ixRMmax2,jpe,2)
   if(((p->mpiupperb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpisendmod(nvar,var3,ixrmmin,ixrmmax,p->jpe,1,2,p);
   //! Wait for messages to arrive
   //call MPI_WAITALL(nmpirequest,mpirequests,mpistatus,ierrmpi)
   /*request.Waitall(gnmpirequest,gmpirequest);*/

   //! Copy buffer received from right (2) physical cells into left ghost cells
   //if(mpilowerB(1) .or. periodic)call mpibuffer2var(2,nvar,var,ixLGmin1,&
   //   ixLGmin2,ixLGmax1,ixLGmax2)
   if(((p->mpilowerb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpibuffer2varmod(1,nvar,var3,ixlgmin,ixlgmax,2,p);
   //! Copy buffer received from left (1) physical cells into right ghost cells
   //if(mpiupperB(1) .or. periodic)call mpibuffer2var(1,nvar,var,ixRGmin1,&
   //   ixRGmin2,ixRGmax1,ixRGmax2)
   if(((p->mpiupperb[2])==1) ||  ((p->boundtype[0][2][0])==0)||  ((p->boundtype[0][0][0])==2))
             mpibuffer2varmod(0,nvar,var3,ixrgmin,ixrgmax,2,p);
}


#endif


}




//!=============================================================================
//subroutine mpireduce(a,mpifunc)
//
//! reduce input for one PE 0 using mpifunc
//
//include 'mpif.h'
//
//double precision :: a, alocal
//integer          :: mpifunc, ierrmpi
//!----------------------------------------------------------------------------
//alocal = a
//call MPI_REDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,0,MPI_COMM_WORLD,&
//   ierrmpi)
void mpireduce(real *a, void *mpifunc)
{
  real alocal;

   alocal=*a;
   /*comm.Reduce(&alocal, a, 1, MPI_DOUBLE_PRECISION, mpifunc, 0);*/

}



//!==============================================================================
//subroutine mpiallreduce(a,mpifunc)
//
//! reduce input onto all PE-s using mpifunc
//
//include 'mpif.h'
//
//double precision :: a, alocal
//integer          :: mpifunc, ierrmpi
//!-----------------------------------------------------------------------------
//alocal = a
//call MPI_ALLREDUCE(alocal,a,1,MPI_DOUBLE_PRECISION,mpifunc,MPI_COMM_WORLD,&
//   ierrmpi)
void mpiallreduce(real *a, void *mpifunc)
{
   real alocal;

   alocal=*a;
   /*comm.Allreduce(&alocal, a, 1, MPI_DOUBLE_PRECISION, mpifunc);*/

}

//update viscosity term
void mpivisc( int idim,struct params *p, real *var1, real *var2, real *var3)
{
   //comm.Barrier();
   int i,n;
   int i1,i2,i3;
   int bound;
   i3=0;
   #ifdef USE_SAC3D
   n=(p->n[0])*(p->n[1])*(p->n[2]);
   switch(idim)
   {
               case 0:
                    n/=(p->n[0]);
                    break;
               case 1:
                    n/=(p->n[1]);
                    break;
               case 2:
                    n/=(p->n[2]);
                    break;               
               }
   #else
   n=(p->n[0])*(p->n[1]);

   switch(idim)
   {
               case 0:
                    n/=(p->n[0]);
                    break;
               case 1:
                    n/=(p->n[1]);
                    break;           
               }   
   #endif
   
   n*=2; //multiply up for all four boundaries
   switch(idim)
   {
   case 0:

     mgpuneighbours(0, p);

     gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
     
        if((p->mpiupperb[idim])==1  ) gnmpirequest++;
        if((p->mpiupperb[idim])==1 )
{ 
//printf("vis1bupper proc %d %d %d  %d %d %d %d\n",p->ipe,p->jpe,p->phpe,p->mpiupperb[idim],p->mpilowerb[idim],n,100*(p->jpe)+10*(idim+1));
//gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferr[0],n,MPI_DOUBLE_PRECISION,p->jpe,100*(p->jpe)+10*(idim+1));
;//gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferr[0],n,MPI_DOUBLE_PRECISION,p->pjpe[idim],MPI_ANY_TAG);

}
//comm.Barrier();
        if((p->mpilowerb[idim])==1  ) gnmpirequest++;
//comm.Barrier();
        if((p->mpilowerb[idim])==1  )
        {
//printf("vis1blower proc %d %d %d  %d %d %d %d\n",p->ipe,p->pjpe[idim],p->phpe[idim],p->mpiupperb[idim],p->mpilowerb[idim],n,100*((p->hpe))+10*(idim+1)+1);
 //        gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferl[0],n,MPI_DOUBLE_PRECISION,p->hpe,100*(p->hpe)+10*(idim+1)+1);
 //        gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferl[0],n,MPI_DOUBLE_PRECISION,p->phpe[idim],MPI_ANY_TAG);

        }
        /*comm.Barrier();*/
          
        if((p->mpiupperb[idim])==1  ) 
{

//printf("vis1a upperproc %d %d %d  %d %d   %d %d\n",p->ipe,p->jpe,p->hpe,p->mpiupperb[idim],p->mpilowerb[idim],n, 100*(p->ipe)+10*(idim+1)+1);  
/*comm.Rsend(gmpisrcbufferr[0], n, MPI_DOUBLE_PRECISION, p->jpe,100*(p->ipe)+10*(idim+1)+1 );*/
;//comm.Rsend(gmpisrcbufferr[0], n, MPI_DOUBLE_PRECISION, p->pjpe[idim], MPI_ANY_TAG);


}
          
//comm.Barrier();
        if((p->mpilowerb[idim])==1  )
{
//  printf("vis1a lowerproc %d %d %d  %d %d   %d %d\n",p->ipe,p->jpe,p->hpe,p->mpiupperb[idim],p->mpilowerb[idim],n,  100*(p->ipe)+10*(idim+1));
 /*comm.Rsend(gmpisrcbufferl[0], n, MPI_DOUBLE_PRECISION, p->hpe,  100*(p->ipe)+10*(idim+1));*/
;//comm.Rsend(gmpisrcbufferl[0], n, MPI_DOUBLE_PRECISION, p->phpe[idim], MPI_ANY_TAG);


}

 //printf("waiting %d\n",p->ipe);

 //    comm.Barrier();
 //      printf("waiting AFTERB %d\n",p->ipe);
 //       request.Waitall(gnmpirequest,gmpirequest);

//comm.Barrier();
  
        //copy data from buffer to the viscosity data in temp2
        //organise buffers so that pointers are swapped instead

        #ifdef USE_SAC3D
   //tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tgtbufferR1(1,ixFlo2:ixFhi2,&
   //   ixFlo3:ixFhi3) !right, upper R
   //tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2,ixFlo3:ixFhi3)=tgtbufferL1(1,ixFlo2:ixFhi2,&
   //   ixFlo3:ixFhi3) !left, lower  L
         for(i2=1;i2<((p->n[1])+2);i2++ )
                  for(i3=1;i3<((p->n[2])+2);i3++ )
         {
          //i1=(p->n[0])+1;
          //bound=i1;
          //var1[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[0][i2+i3*((p->n[1])+2)];
          //var1[encode3p2_sacmpi (p,0, i2, i3, tmpnui)]=gmpitgtbufferl[0][i2+i3*((p->n[1])+2)];
          i1=0;
          for(bound=0; bound<2; bound++)
            var1[sacencodempivisc0(p,i1,i2,i3,bound+2,idim)]=gmpitgtbufferr[0][i2+i3*((p->n[1])+2)+bound*((p->n[1])+2)*((p->n[2])+2)];
          for(bound=0; bound<2; bound++)
             var1[sacencodempivisc0(p,i1,i2,i3,bound,idim)]=gmpitgtbufferl[0][i2+i3*((p->n[1])+2)+bound*((p->n[1])+2)*((p->n[2])+2)];


          }




        #else
        //tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2)=tgtbufferR1(1,ixFlo2:ixFhi2) !right, upper R
        //tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2)=tgtbufferL1(1,ixFlo2:ixFhi2) !left, lower  L
                  for(i2=1;i2<((p->n[1])+2);i2++ )
         {
          //i1=(p->n[0])+1;
         
         // var1[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[0][i2];
          //var1[encode3p2_sacmpi (p,0, i2, i3, tmpnui)]=gmpitgtbufferl[0][i2];

          i1=0;
          for(bound=0; bound<2; bound++)
            var1[sacencodempivisc0(p,i1,i2,i3,bound+2,idim)]=gmpitgtbufferr[0][i2+bound*((p->n[1])+2)];
          for(bound=0; bound<2; bound++)
             var1[sacencodempivisc0(p,i1,i2,i3,bound,idim)]=gmpitgtbufferl[0][i2+bound*((p->n[1])+2)];


          }

 
         
        #endif

        

     break;
     
        case 1:
     mgpuneighbours(1, p);

     gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
     




        if((p->mpiupperb[idim])==1  ) gnmpirequest++;
        if((p->mpiupperb[idim])==1 )
;//gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferr[1],n,MPI_DOUBLE_PRECISION,p->jpe,100*(p->jpe)+10*(idim+1));
        if((p->mpilowerb[idim])==1  ) gnmpirequest++;
        if((p->mpilowerb[idim])==1  )
        ;// gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferl[1],n,MPI_DOUBLE_PRECISION,p->hpe,100*(p->hpe)+10*(idim+1)+1);
        /*comm.Barrier();*/
          
        if((p->mpiupperb[idim])==1  ) 
/*comm.Rsend(gmpisrcbufferr[1], n, MPI_DOUBLE_PRECISION, p->jpe,100*(p->ipe)+10*(idim+1)+1 );*/
          
        if((p->mpilowerb[idim])==1  )
 /*comm.Rsend(gmpisrcbufferl[1], n, MPI_DOUBLE_PRECISION, p->hpe,  100*(p->ipe)+10*(idim+1));*/

     /*comm.Barrier();
 
        request.Waitall(gnmpirequest,gmpirequest);*/














      #ifdef USE_SAC3D
  //tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1,ixFlo3:ixFhi3)=tgtbufferR2(ixFlo1:ixFhi1,1,&
  //    ixFlo3:ixFhi3) !right, upper R
  // tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1,ixFlo3:ixFhi3)=tgtbufferL2(ixFlo1:ixFhi1,1,&
   //   ixFlo3:ixFhi3) !left, lower  L
         for(i1=1;i1<((p->n[0])+2);i1++ )
                  for(i3=1;i3<((p->n[2])+2);i3++ )
         {
          //i2=(p->n[1])+1;
         
          //var2[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[1][i1+i3*((p->n[0])+2)];
          //var2[encode3p2_sacmpi (p,i1, 0, i3, tmpnui)]=gmpitgtbufferl[1][i1+i3*((p->n[0])+2)];

          i2=0;
          for(bound=0; bound<2; bound++)
            var2[sacencodempivisc1(p,i1,i2,i3,bound+2,idim)]=gmpitgtbufferr[0][i1+i3*((p->n[0])+2)+bound*((p->n[0])+2)*((p->n[3])+2)];
          for(bound=0; bound<2; bound++)
             var2[sacencodempivisc1(p,i1,i2,i3,bound,idim)]=gmpitgtbufferl[0][i1+i3*((p->n[0])+2)+bound*((p->n[0])+2)*((p->n[3])+2)];




          }


     #else
       // tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1)=tgtbufferR2(ixFlo1:ixFhi1,1) !right, upper R
   //tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1)=tgtbufferL2(ixFlo1:ixFhi1,1) !left, lower  L
                  for(i1=1;i1<((p->n[0])+2);i1++ )
         {
          //i2=(p->n[1])+1;
         
          //var2[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[1][i1];
          //var2[encode3p2_sacmpi (p,i1, 0, i3, tmpnui)]=gmpitgtbufferl[1][i1];


          i2=0;
          for(bound=0; bound<2; bound++)
            var2[sacencodempivisc1(p,i1,i2,i3,bound+2,idim)]=gmpitgtbufferr[0][i1+bound*((p->n[0])+2)];
          for(bound=0; bound<2; bound++)
             var2[sacencodempivisc1(p,i1,i2,i3,bound,idim)]=gmpitgtbufferl[0][i1+bound*((p->n[0])+2)];




          }

      #endif
     break;
     #ifdef USE_SAC3D
        case 2:
      mgpuneighbours(2, p);

     gnmpirequest=0;  
  for(i=0; i<2; i++)
     gmpirequest[i]=MPI_REQUEST_NULL;
   


        if((p->mpiupperb[idim])==1  ) gnmpirequest++;
       /* if((p->mpiupperb[idim])==1 )
gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferr[2],n,MPI_DOUBLE_PRECISION,p->jpe,100*(p->jpe)+10*(idim+1));*/
        if((p->mpilowerb[idim])==1  ) gnmpirequest++;
       /* if((p->mpilowerb[idim])==1  )
         gmpirequest[gnmpirequest]=comm.Irecv(gmpitgtbufferl[2],n,MPI_DOUBLE_PRECISION,p->hpe,100*(p->hpe)+10*(idim+1)+1);
        comm.Barrier();*/
          
      /*  if((p->mpiupperb[idim])==1  ) 
comm.Rsend(gmpisrcbufferr[2], n, MPI_DOUBLE_PRECISION, p->jpe,100*(p->ipe)+10*(idim+1)+1 );
          
        if((p->mpilowerb[idim])==1  )
 comm.Rsend(gmpisrcbufferl[2], n, MPI_DOUBLE_PRECISION, p->hpe,  100*(p->ipe)+10*(idim+1));

     comm.Barrier();
 
        request.Waitall(gnmpirequest,gmpirequest);  */



   //  tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFhi3+1)=tgtbufferR3(ixFlo1:ixFhi1,&
   //   ixFlo2:ixFhi2,1) !right, upper R

   //tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2,ixFlo3-1)=tgtbufferL3(ixFlo1:ixFhi1,&
   //   ixFlo2:ixFhi2,1) !left, lower  L
        for(i1=1;i1<((p->n[0])+2);i1++ )
                  for(i2=1;i2<((p->n[1])+2);i2++ )
         {
         // i3=(p->n[2])+1;
         
        //  var3[encode3p2_sacmpi (p,i1, i2, i3, tmpnui)]=gmpitgtbufferr[2][i1+i2*((p->n[0])+2)];
         // var3[encode3p2_sacmpi (p,i1, i2, 0, tmpnui)]=gmpitgtbufferl[2][i1+i2*((p->n[0])+2)];


          i3=0;
          for(bound=0; bound<2; bound++)
            var3[sacencodempivisc2(p,i1,i2,i3,bound+2,idim)]=gmpitgtbufferr[0][i1+i2*((p->n[0])+2)+bound*((p->n[0])+2)*((p->n[1])+2)];
          for(bound=0; bound<2; bound++)
             var3[sacencodempivisc2(p,i1,i2,i3,bound,idim)]=gmpitgtbufferl[0][i1+i2*((p->n[0])+2)+bound*((p->n[0])+2)*((p->n[1])+2)];




          }

     
     break;
     #endif
   
  }
  // comm.Barrier();
}










#endif




