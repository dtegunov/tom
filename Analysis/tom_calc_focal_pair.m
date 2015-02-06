function Dz2=tom_calc_focal_pair(Dz1,pix_size,voltage,pixs,Cs,df_range)


%calc 1 ctf

if (nargin < 6)
    df_range=[-3.2 -0.05 -8];
end;

voltagest=(voltage.* 1000)*(1+(voltage.*1000)/1022000);

Dz1_m=Dz1*10^(-6);
ctf1=tom_ctf(Dz1,pix_size,voltage,pixs,Cs);
lambda=sqrt(150.4/voltagest)*10^-10;


ctf_zero = sqrt(lambda*abs(Dz1_m)*10^18);

zz=1;
for i=df_range(1):df_range(2):df_range(3)
    ctf2=tom_ctf(i,pix_size,voltage,pixs,Cs);
    d_fun=diff(ctf2);
    [index b]=tom_crossing(d_fun,[],0);
    res(zz)=((pixs./2)./b(2)).*(pix_size.*2);
    dz_samp(zz)=i;
    zz=zz+1;
end;

[pointidx, pointcoords, distance] = tom_nearestpoint(ctf_zero,res);

Dz2=dz_samp(pointidx);
 
disp(' ');