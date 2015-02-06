function [all_fsc,res_vect,sum_vect]=tom_av2_xmipp_all_fsc(basename_fold,basename_res,extend_res,crit,coef_thr)
%TOM_AV2_XMIPP_ALL_FSC reads all fsc according 2 given wildcard
%
%   all_fsc=tom_av2_xmipp_all_fsc(wildcard)
%
%PARAMETERS
%
%  INPUT
%   basename_fold          basename of fold    
%   basename_res           basename of resolution file  
%   extend_res             extensoin of resolution file  
%   crit                   (opt.) criteria for res mesurement (e.g. 0.5 0.3 0.14)                     
%   coef_thr               (-1)   threshold for setting fsc 2 0
%
%  OUTPUT
%   all_fsc           all fsc as 2d array
%   res_vect          resolution for given crit 
%   sum_vect          sum of fsc coeffs
%EXAMPLE
%
% all_fsc=tom_av2_xmipp_all_fsc('ProjMatch/run_org/Iter_','Iter_','_resolution.fsc');
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_angular_plot
%
%   created byFB 10/31/12
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


if (nargin<4)
    crit='';
end;
if (nargin<5)
    coef_thr=-1;
end;
res_vect=[];

d=dir([basename_fold '*']);
if (isempty(d))
    disp('no fsc files found! ...check inputs')
    all_fsc='';
    return;
end;

for i=1:length(d)
    name=[basename_fold num2str(i) '/' basename_res num2str(i) extend_res];
    try
        dat=importdata(name);
    catch Me
        disp(' ');
        disp([name ' not found stop']);
        return;
    end;
    tmp_dat=dat.data(:,2);
    idx=find(tmp_dat<coef_thr);
    if (isempty(idx)==0)
        tmp_dat(idx(1):end)=0;
        dat.data(idx(1):end,2)=0;
    end
    all_fsc(:,i)=tmp_dat;
    
    if (isempty(crit)==0)
        [res_vect(i) sum_vect(i)]=calc_res(dat.data,crit); %.data,dat.data(end,4),size(dat.data,1),crit);
    end;
    disp(name);
    
    
end;


function [v,csum]=calc_res(data,crit)


% for i=1:num_of_sh
%     res=(num_of_sh./i).*ny;
%     data(i,1)=1./res;
%     if (isnan(mat(i,9)) )
%         data(i,2)=0.999;
%     else
%         data(i,2)=mat(i,9);
%     end
%     data(i,3)=-1;
%     data(i,4)=res;
% end;


v=data(1,4); %ini with lowest res
for i=1:size(data(:,2),1)
    if data(i,2)<crit
        if (i>1)
            v1=data(i-1,4);
        else
            v1=data(i,4);
        end;
        v2=data(i,4);
        v=(v1+v2)./2;
        break;
    end;
end;

idx=find(data(:,2)>0);
csum=sum(data(idx,2));









