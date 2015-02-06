function [st_doc mat]=tom_av2_xmipp_analyse_cc_ctf(doc_filename,cftDat_file,f_doc_struct_ext)



%read lists
fprintf('\n%s ', ['Reading  ' doc_filename ':']);  
st_doc=tom_xmippdocread(doc_filename);

fprintf('%s \n', ['...done!']); 


mat=zeros(6,length(st_doc));
zz_not=0;
zz_p=1;

five_p=length(st_doc)./20;
tic;

for i=1:length(st_doc)

    [a b]=unix(['cat ' cftDat_file ' | grep ' st_doc(i).name]);
    [c d]=strtok(b,' ');
    d=deblank(d); d=strtrim(d);
    ctf_param_name=d;
    
   
    try
       ctf_dat=importdata(ctf_param_name);
       
       %transfer data
       st_doc(i).sampling_rate=ctf_dat.data(1);
       st_doc(i).voltage=ctf_dat.data(2);
       st_doc(i).defocusU=ctf_dat.data(3).*1e-4;
       st_doc(i).defocusV=ctf_dat.data(4).*1e-4;
       st_doc(i).azimuthal_angle=ctf_dat.data(5);
       st_doc(i).spherical_aberration=ctf_dat.data(6);
       st_doc(i).chromatic_aberration=ctf_dat.data(7);
       
       mat(1,i)=st_doc(i).defocusU;
       mat(2,i)=st_doc(i).defocusV;
       mat(3,i)=st_doc(i).azimuthal_angle;
       mat(4,i)=abs(st_doc(i).defocusU-st_doc(i).defocusV);
       mat(5,i)=mean([st_doc(i).defocusU st_doc(i).defocusV]);
       mat(6,i)=st_doc(i).pmax_sump;
        
   
   catch
       disp([ctf_param_name ' not readable !']);
       st_doc(i).sampling_rate=-1;
       st_doc(i).voltage=-1;
       st_doc(i).defocusU=-1;
       st_doc(i).defocusV=-1;
       st_doc(i).azimuthal_angle=-1;
       st_doc(i).spherical_aberration=-1;
       st_doc(i).chromatic_aberration=-1;
       
       mat(1,i)=-1;
       mat(2,i)=-1;
       mat(3,i)=-1;
       mat(4,i)=-1;
       mat(5,i)=-1;
       mat(6,i)=-1;
       zz_not=zz_not+1;
       
   end;
  
   if (mod(i,(five_p.*zz_p))==0)
        toc;
        disp([num2str(zz_p.*5) ' % done ']);
        tic;
        zz_p=zz_p+1;
   end;
   
   if (mod(i,100)==0)  
       fprintf('%s','.');  
   end;
   
end;

fprintf('%s \n', [ 'done!']); 
disp([num2str(zz_not) ' of ' num2str(length(st_doc)) ' ctf-param files not found ']);


if (nargin > 3)
    save(f_doc_struct_ext,'st_doc');
end;
    


% [b all]=unix(['awk ''{split($2,a,"/");  for (i=2;i<length(a);i++) { printf "/%s" ,a[i]  }  split(a[length(a)],d,"_"); printf "/%s; \n" ,d[1]  }'' ' cftDat_file]);
% 
% C = textscan(all,'%s','Delimiter',';');
% C=C{1};
% u_path=unique(C);

%grep "defocusU\|defocusV"
%/fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/090721_p47f11/log/rec/parts_low_128_ctfDat/ctf_*.ctfparam;


