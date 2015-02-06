function [stack rots shifts]=tom_av2_create_alg_phantom(sz,number)


tmpl=tom_rectanglemask(zeros(sz),[(sz(1)./3) (sz(2)./3)],1,[(sz./2) (sz./2) 1]);
%tmpl=tmpl==0;
mask=tom_spheremask(ones(sz),sz(1)./8,0,[((sz(1)./2)+(sz(1)./6)) ((sz(1)./2)+(sz(1)./6)) 1]);
%maks=mask==0;
tmpl=double((tmpl+mask)>=0.5);
tom_emwrite('p_ref',tmpl);
%rot_lapaloma=[0 48 93 142 95 5 45 90 135 90 0 45 90 135 90 0 45 90 135 90];

shifts =[   24.4772   21.0542
   11.9758    9.7994
   22.6159    5.3877
   14.9193   23.1131
   11.6001   23.3865
   21.8029    0.2436
   30.2475   19.6247
    0.1009   25.5027
    5.7113   16.9408
   17.5377    1.8636
   13.3153    5.9662
    2.3930    9.9212
   31.3833   17.7637
   22.1299    7.7331
   29.9044    4.1217
    9.5118   20.7111
   29.5293    7.7351
   10.6332   15.2561
   22.5897    7.6784
   27.6863   13.1517];


rot_lapaloma=[0 90 90 0 0 0 90 90 0 0 0 90 90 0 0 0 90 90 0 0];

rots=[  259.5294
  291.7471
  133.4519
  292.9396
  358.1528
  233.8803
  235.4937
  281.8541
  231.0540
   78.7477
  211.5314
   23.0079
  339.8706
  355.8692
  291.5329
  247.2577
  166.9532
  237.6563
  168.7584
  258.1929
];

for i=1:number
    rot=(rand(1).*360);
    %rot=rot_lapaloma(i);
%    rots(i)=rot;
    %rot=0;
   % rot=rot_lapaloma(i);
%    shift=[rand(1).*(sz(1)./5) rand(1).*(sz(2)./5)];
%    shift=[rand(1).*(sz(1)./5) rand(1).*(sz(2)./5)];
    %shift=[0 0];
%    shifts(i,:)=shift;
    out=tom_rotate(tmpl,rots(i));
    out=tom_shift(out,shifts(i,:));
    
    %out=out+rand(1).*100;
    stack(:,:,i)=out;
end;


%and back!
nr=zeros(sz);
for i=1:number
    im=stack(:,:,i);
    out=tom_rotate(im,-rots(i));
    nr=nr+out;
    
end;

disp('end');