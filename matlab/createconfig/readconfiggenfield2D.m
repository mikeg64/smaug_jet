
%read a configuration file
%getpicttest  3D version
% Read the npict-th picture from 1 or more files
filename='/data/cs1mkg/smaug_spicule1/spicule4b1_3_3d/zerospic1__711000.out';

   fid=fopen(trim(filename));
   %fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   headline=trim(setstr(fread(fid,79,'char')'));
   it=fread(fid,1,'integer*4'); time=fread(fid,1,'float64');
 
   ndim=fread(fid,1,'integer*4');
   neqpar=fread(fid,1,'integer*4'); 
   nw=fread(fid,1,'integer*4');
   nx=fread(fid,2,'integer*4');
   
   nxs=nx(1)*nx(2);%*nx(3);
   varbuf=fread(fid,7,'float64');
   
   gamma=varbuf(1);
   eta=varbuf(2);
   g(1)=varbuf(3);
   g(2)=varbuf(4);
   %g(3)=varbuf(5);
   
   
   varnames=trim(setstr(fread(fid,79,'char')'));
   
   for idim=1:ndim
      X(:,idim)=fread(fid,nxs,'float64');
   end
   
   for iw=1:nw
      %fread(fid,4);
      w(:,iw)=fread(fid,nxs,'float64');
      %fread(fid,4);
   end
   
   nx1=nx(1);
   nx2=nx(2);
   %nx3=nx(3);
   
   xx=reshape(X(:,1),nx1,nx2);
   yy=reshape(X(:,2),nx1,nx2);
   %zz=reshape(X(:,3),nx1,nx2,nx3);
   
   
 
  % extract variables from w into variables named after the strings in wnames
wd=zeros(nw,nx1,nx2);
for iw=1:nw
  
     tmp=reshape(w(:,iw),nx1,nx2);
     wd(iw,:,:)=tmp;
end


%w=tmp(iw);
  

clear tmp; 
   
   
   fclose(fid);
    





    



newfilename='2D_2048_1024_8_8_asc.ini';

% 4Mmx12.5Mm
%extrapolating atmospheric model to greater than 6Mm
%sets the correct minimum and maximum so that the ghost cells are also
%included
simparams=sim_params;
simgridinfo=sim_gridinfo;
simdata=sim_data;

simdata.w=zeros(nx1,nx2,nw);


for(i=0:10)
    simdata.w(:,:,i)=wd(i,:,:);
end





simparams.unique_identifier="2dhydro 2048x1024";

consts.mu=0.6e0; %magnetic permeability
consts.R=8.31e3
consts.fgamma=1.66666667e0
consts.ggg=274.0e0 % acceleration due to gravity on the sun
consts.mu=4*pi/1.0e7

  it=0;
       time=0;
       ndim=2;
       neqpar=6 ;%neqpar=7;
       nw=10; %nw=13;
       nx1=2048;
       nx2=1024;
       nx3=1;
       
       ng1=0;
       ng2=0;
       ng3=0;


      gamma=1.666667;
       eta=0;
       g2=-274;
       g1=0;
       g3=0;
       
       %xmin=133333.33;
       ymin=199219.0;
       %ymin=1953.1;
       xmin=39687.5;
       %zmin=1953.1;
       zmin=39687.5;
       %xmax=5955555.6e0;
       %ymax=12.8496e6;
       ymax=8.0e6;
       xmax=8.0e6;
       zmax=8.0e6;
       %xmax=2.559984e6;
       %zmax=2.559984e6;
       
       dx=(xmax-xmin)/(nx1-1-2*ng1);
       dy=(ymax-ymin)/(nx2-1-2*ng2);
       %dz=(zmax-zmin)/(nx3-1-2*ng3);
       
       xmin=xmin-2*dx;
       xmax=xmax+2*dx;
       ymin=ymin-2*dy;
       ymax=ymax+2*dy;
       

       dx=(xmax-xmin)/(nx1-1-2*ng1);
       dy=(ymax-ymin)/(nx2-1-2*ng2);
       %dz=(zmax-zmin)/(nx3-1-2*ng3);      
       
%        xx=zeros(nx1,nx2);
%        yy=zeros(nx1,nx2);
%        zz=zeros(nx1,nx2);
%        
%        
%        for i=1:nx1
%            for j=1:nx2
%                %for k=1:nx3
%                rheight(j)=ymin-(ng2*dy)+dy*(j-1);
%                xx(i,j)=xmin-(ng1*dx)+dx*(i-1);
%                yy(i,j)=ymin-(ng2*dy)+dy*(j-1);
%                %zz(i,j,k)=zmin-(ng3*dz)+dz*(k-1);
%                %end
%            end
%        end
%        
       
       




        simparams.current_iteration=it;
        simparams.current_time=time;
        simparams.dimensionality=ndim;
        simparams.domain_dimensions=[nx1;nx2;nx3];
%         simparams.domain_left_edge=[0;0; 0.0];
%         simparams.domain_right_edge=[0; 0; 0];
        simparams.eta=eta;
%        field_ordering=1;
        simparams.gamma=gamma;
        simparams.gravity0=g1;
        simparams.gravity1=g2;
        simparams.gravity2=g3;
        
        
        
        
        simparams.domain_left_edge(1)=xx(1,1);
        simparams.domain_left_edge(2)=yy(1,1);
        simparams.domain_left_edge(3)=0;
        
        simparams.domain_right_edge(1)=xx(nx1,1);
        simparams.domain_right_edge(2)=yy(1,nx2);
        simparams.domain_right_edge(3)=0;

        
        simgridinfo.grid_dimensions(1)=nx1;
        simgridinfo.grid_dimensions(2)=nx2;
        simgridinfo.grid_dimensions(3)=nx3;
        
        simgridinfo.ndimensions=2;
        









[simparams, simgridinfo, simdata]=generatefield_verttube(simparams, simgridinfo, simdata, 'vertfieldtube');         
      
      
      
 writesac2D(newfilename, simparams, simgridinfo, simdata, 'ascii');
        