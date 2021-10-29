function [simparams, simgridinfo, simdata]=generatefield_verttube(simparams, simgridinfo, simdata, mode)

%Generate magnetic field configuration
%simple fluxtube using self similarity and hydrostatic pressure correction



%Generate field
% calculate hydrostatic pressure update


%cases
%1.  Vertical field
%2. Horizontal field
%3. Inclined (off vertical field)
%4. Flux tube array

if strcmp(mode,'vertfieldtube')
    
    
    %tube width =100km
    %footpoint intensity 1kG ala thin photospheric flux tubes
    nb=1;
    mb=1;
    
    p0=5226.1d0; %read the pressure from the valIIc data file corresponding to minimum height of the mode
    
    %swap y and x directions
    nx3=simgridinfo.grid_dimensions(1);
    nx1=simgridinfo.grid_dimensions(2);
    
    dz=-(simparams.domain_left_edge(1)-simparams.domain_right_edge(1))/simgridinfo.grid_dimensions(1);
    dx=-(simparams.domain_left_edge(2)-simparams.domain_right_edge(2))/simgridinfo.grid_dimensions(2);

    bx=zeros(nx1,nx3);
    bz=zeros(nx1,nx3);
    b0z=zeros(nx1,1);
    
    xf=zeros(nx1,nx3);
    
Bmax=0.15  ; %mag field Tesla
%Bmax=0.0003;
%Bmin=0.0006d0  ; %mag field Tesla
Bmin=0.0002  ; %mag field Tesla

 b0zz=Bmax;   
    d_z=1.5; % width of Gaussian in Mm
    z_shift= 0.0e6; % shift in Mm
    %archshift=1.0e6;
    archshift=0.0;
    archlength=8.0e6;
    A=0.45; % amplitude
    scale=1.0e6;
    barch0=0.0;  %0.05;
    b0z_top=0.08;
    l=4e6; %length of arcade
   
    f0=2.0e6; %tube opening factor

    Ab0z=20.e0; % bz - amplitude
    mu=4.0*pi*1.0e-7
 

gamma=1.66666667;
ggg=-274.0;
   

xr=0.1e6;
yr=0.1e6;

R2=(xr.^2+yr.^2);

A=R2/2;


  
    x=zeros(nx1,1);
    z=zeros(nx3,1);
    

    for k=1:nx3
        z(k)=simparams.domain_left_edge(2)+dz*(k-1);
    end
    
    for i=1:nx1
        x(i)=simparams.domain_left_edge(1)+dx*(i-1);
        %b0z(i)=(par4((x(i)/scale-z_shift),d_z,A)).^2;
        %b0z(i)=(par3((x(i)/scale-z_shift),d_z,A));
        b0z(i)=par4((x(i)/scale-z_shift),d_z,A);
    end
    
    dzb=(simparams.domain_right_edge(1)-simparams.domain_left_edge(1))/nb;
    %dzb=(simparams.domain_right_edge(3)-simparams.domain_left_edge(3))/mb;
    
    
    
	bnmin=min(b0z);
	bnmax=max(b0z);
    
    
    for i=1:nx1
	b0z(i)=((Bmax-Bmin)./(bnmax-bnmin)).*(b0z(i)-bnmin)+Bmin;
    end
    
    %b0z=b0z./max(b0z);  
    %b0z=Ab0z.*b0z+b0z_top;
    
    dbz=deriv1(b0z,x,1);


	%y=y-(max(y)-min(y))/2 ;%+5000.d0
	%z=z-(max(z)-min(z))/2 ;%+5000.d0

  xold=x;
  %yold=y;
  zold=z;
  %z=zold;  
    
 ib=0;
 jb=0;
  %  for ib=0:nb-1
  %      for jb=0:mb-1





            
          %  if nb>1
          %      ybp=simparams.domain_left_edge(2)+(ib+1)*nx2*dy*(simgridinfo.grid_dimensions(2)/((nb+1)))-nx2*dy*(simgridinfo.grid_dimensions(2)/(2*(nb+1)));
          %  else
          %      ybp=simparams.domain_left_edge(2)+(ib+1)*nx2*dy*(1/((nb+1)));                
          %  end
            
          %  if mb>1
          %      zbp=simparams.domain_left_edge(3)+(jb+1)*nx3*dz*(simgridinfo.grid_dimensions(3)/((mb+1)))-nx3*dz*(simgridinfo.grid_dimensions(3)/(2*(mb+1)));                
          %  else
          %      zbp=simparams.domain_left_edge(3)+(jb+1)*nx3*dz*(1/((mb+1)));
          %  end
            
  
            z=zold-(max(z)+(dzb/nx3))/2.d0;
          %  y=yold-max(y)/2.d0;

            
           % z=zold-((jb-1)*dzb/2)-dzb;
           %z=zold+2*d_z-(dzb)/2;
            %y=yold-((ib-1)*dyb/2)-dyb;           
            
            
            
            
for k=1:nx3
for i=1:nx1

%f=b0z(i)*sqrt((y(j)).^2+(z(k)).^2);

%xf(i,j,k)=(par4(f,f0,0.5)).^2;


f=(z(k).^2)./R2;
xf(i,k)=exp(-f);





end
end
xf=xf./max(max(xf));
%b0zz=0.15;

for k=1:nx3 
for i=1:nx1

% bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
% bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j)-ybp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k));
% by(i,j,k)=by(i,j,k)-dbz(i)*(y(k)-zbp)/sqrt((x(j)-ybp).^2+(y(k)-zbp).^2)*xf(i,j,k);
%bz(i,j,k)=bz(i,j,k)+(b0z(i)/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
%bx(i,j,k)=bx(i,j,k)-(dbz(i)*(x(j))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k));
%by(i,j,k)=by(i,j,k)-dbz(i)*(y(k))/sqrt((x(j)).^2+(y(k)).^2)*xf(i,j,k);


bz(i,k)=bz(i,k)+b0zz.*xf(i,k);
bx(i,k)=bx(i,k);
%by(i,j,k)=by(i,j,k);

l=archlength;
%expressions for arcade
 bz(i,k)=bz(i,k)+barch0*cos(pi*(x(i)-archshift)/(l))*exp(-pi*z(k)/(l));
 bx(i,k)=bx(i,k)+barch0*sin(pi*(x(i)-archshift)/(l))*exp(-pi*z(k)/(l));



end
end  

disp(b0zz);

zold=z;
%yold=y;


                       
       %  end  %end of loop
 %    end




%     
 %end   %building fluxtube

bb=(bx.^2+bz.^2)/2/mu;


% ************** convert to VAC magnetic field
bz=bz./sqrt(mu);
bx=bx./sqrt(mu);
%by=by./sqrt(mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%checked to here


%compute magnetostatics pressure correction
dbzdz=zeros(nx1,nx3);
dbxdz=zeros(nx1,nx3);
%dbydz=zeros(nx1,nx2,nx3);

dbzdx=zeros(nx1,nx3);
dbxdx=zeros(nx1,nx3);
%dbydx=zeros(nx1,nx2,nx3);

%dbzdy=zeros(nx1,nx2,nx3);
%dbxdy=zeros(nx1,nx2,nx3);
%dbydy=zeros(nx1,nx2,nx3);


Bvarix=zeros(nx1,nx3);
%Bvariy=zeros(nx1,nx2,nx3);
Bvari=zeros(nx1,nx3);

dpdz=zeros(nx1,nx3);

rho1=zeros(nx1,nx3);

Bvaridz=zeros(nx1,nx3);
Bvaridz1=zeros(nx1,nx3);
dBvardz=zeros(nx1,nx3);
br=zeros(nx1,nx3);

%bxby=zeros(nx1,nx2,nx3);
%dbxbydy=zeros(nx1,nx2,nx3);
bxbz=zeros(nx1,nx3);
bxdbzdx=zeros(nx1,nx3);
bxdbzdz=zeros(nx1,nx3);
dbxbzdz=zeros(nx1,nx3);
dbxbzdx=zeros(nx1,nx3);
%bxby=zeros(nx1,nx2,nx3);
%bybz=zeros(nx1,nx2,nx3);
%dbybzdy=zeros(nx1,nx2,nx3);

for k=1:nx3
%for j=1:nx2

dbxdx(:,k)=deriv1(bx(:,k),x,1);

% dbxdz(:,j,k)=deriv1(bx(:,j,k),x);
% dbydz(:,j,k)=deriv1(by(:,j,k),x);
%end
end



for k=1:nx3
%for j=1:nx2

dbzdx(:,k)=deriv1(bz(:,k),x,1);

% dbxdz(:,j,k)=deriv1(bx(:,j,k),x);
% dbydz(:,j,k)=deriv1(by(:,j,k),x);
%end
end

bxdbzdx=bx.*dbzdx;


%for k=1:nx3
for i=1:nx1
% dbzdx(i,:,k)=deriv1(bz(i,:,k),y);
  dbzdz(i,:)=deriv1(bz(i,:),z,2);
% dbydx(i,:,k)=deriv1(by(i,:,k),y);
%end
end

bxdbzdz=bx.*dbzdz;


%for j=1:nx2
%for i=1:nx1
% dbzdy(i,j,:)=deriv1(bz(i,j,:),z);
% dbxdy(i,j,:)=deriv1(bx(i,j,:),z);
% dbydy(i,j,:)=deriv1(by(i,j,:),z,3); 
%end
%end

divb=dbzdz+dbxdx; %+dbydy;

%bxby=bx.*by;

%%define matlab code here
%check b_field_vertical_tube.pro
%bxby=bx.*by;

%for j=1:nx2
%for i=1:nx1
% dbxbydy(i,j,:)=deriv1(bxby(i,j,:),z,3);
%end
%end


%print,'dBxBydy'


bxbz=bx.*bz;

% %for j=1:nx2
% for i=1:nx1
%  dbxbzdz(i,:)=deriv1(bxbz(i,:),z,2);
% end
% %end
% 
% %for j=1:nx2
% for k=1:nx3
%  dbxbzdx(i,:)=deriv1(bxbz(:,j),x,1);
% end
% %end


%;print,'dBxBzdz'


%bxby=bx.*by;

%for i=1:nx1
%for k=1:nx3
% dbxbydx(i,:,k)=deriv1(bxby(i,:,k),y,2);
%end
%end

%print,'dBxBydx'



%bybz=dblarr(n1,n2,n3)
%dbybzdz=dblarr(n1,n2,n3)
%bybz=by.*bz;

%for j=1:nx2
%for k=1:nx3
% dbybzdz(:,j,k)=deriv1(bybz(:,j,k),x,1);
%end
%end

%************* BEGIN INTEGRATION ****************************
%F=dbxbydy+dbxbzdz;
%G=dbxbydx+dbybzdz;

%F=+dbzbzdz;
F=bx.*dbzdx+bz.*dbzdz;
G=0;

disp('shape F');
disp(size(F));

for i=1:nx1

%for kx=1:nx3
  %for jx=1:nx2
%   sum=inte(reshape(F(i,kx),[i,1]),z(2)-z(1)); 
 sum=inte(reshape(F(i,:),[nx3,1]),z(2)-z(1)); 
   %sum=F(i,kx);
  Bvarix(i)=sum;
% end
end

%for jy=1:nx2
%  for ky=1:nx3
%   sum=inte(reshape(G(i,jy,1:ky),[ky,1]),z(2)-z(1)); 
%  Bvariy(i,jy,ky)=sum;
% end
%end
disp(  i);

%end

%************* END INTEGRATION ****************************

%Bvari=((Bvarix+Bvariy)/2)-(bz.^2)/2;
Bvari=((Bvarix))-(bz.^2)/2;
%Bvari=Bvarix;



%for j=1:nx2
for i=1:nx1
 dpdz(i,:)=deriv1(Bvari(i,:),z,2);
end
%end

bxbz=(bx.*bx+bz.*bz)/2;

for i=1:nx1
%for k=1:nx3
 dbxbzdz(i,:)=deriv1(bxbz(i,:),z,2);
%end
end



% bxbz=bx.*bz;
% 
% for i=1:nx1
% for k=1:nx3
%  dbxbzdx(i,:,k)=deriv1(bxbz(i,:,k),x,2);
% end
% end


% bybz=dblarr(n1,n2,n3)
% dbybzdy=dblarr(n1,n2,n3)
%bybz=by.*bz;

% for i=1:nx1
% for j=1:nx2
%  dbybzdy(i,j,:)=deriv1(bybz(i,j,:),y,3);
% end
% end


%for i=1:nx1
%for j=1:nx2
%for k=1:nx3
% rho1(i,j,k)=(dbxbybzdz(i,j,k)-dbxbzdx(i,j,k)-  dbybzdy(i,j,k)+dpdz(i,j,k))./ggg;
%end
%end
%end

rho1=(dbxbzdz-bxdbzdx-bxdbzdz+dpdz)./ggg;
%rho1=(dbxbzdz-dbxbzdx-  dbybzdy+dpdz)./ggg;

%note the use of transpose operation because we swapped the vertical and x
%directions
rho1=(simdata.w(:,:,1))'+rho1+(simdata.w(:,:,8))';
p=p0+Bvari+(simdata.w(:,:,7))'*(gamma-1.0);


%lower boundary

for ix_1=4:-1:2
  %for ix_2=1:nx2
  for ix_3=1:nx3  
         p_2=rho1(ix_1,ix_3)*ggg;
         p(ix_1-1,ix_3) = (z(2)-z(1))*p_2+p(ix_1,ix_3);
  end  
  %end
 end


%upper boundary

for ix_1=nx1-2:nx1-1
   %for ix_2=1:nx2
   for ix_3=1:nx3   
           p_2=rho1(ix_1,ix_3)*ggg;
           p(ix_1+1,ix_3) = -(z(2)-z(1))*p_2+p(ix_1,ix_3);
   end	   
   %end
end

% minp=min(min(p));
% for k=1:nx3 
% for i=1:nx1
%    if p(i,k)<0
%        p(i,k)=p(i,k)+1.01*abs(p(i,k));
%    end
%     
% end
% end

%update the background energy and magnetic fields
simdata.w(:,:,8)=rho1';
simdata.w(:,:,7)=(p./((gamma-1.0))+0.5*(bx.*bx+bz.*bz))';
simdata.w(:,:,9)=bx';
simdata.w(:,:,10)=bz';


end %if vert tube loop

%rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
%set background density
%set background energy

end








