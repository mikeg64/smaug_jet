xc1=2.0e6;
xc2=800000.0;

qt=1.0;
aa=1.0;

xp=0.0e6:0.1e6:7.0e6;
yp=993465.0;

s_period=100.0;
tdep=1.00;
%tdep=sin(qt*2.0*pi/s_period);

     r1=(xp-xc1).*(xp-xc1);
     r2=(yp-xc2).*(yp-xc2);

 
 s_rad1=500000.0;
 s_rad2=500000.0;



vvv=aa*tdep*exp(-(r1./(s_rad1*s_rad1))-(r2./(s_rad2*s_rad2)));

plot(vvv);