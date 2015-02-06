function [stat error_count]=tom_av2_xmipp_checkstack(in_sel,exp_size,exp_pix_range,norm,subset,verbose)
%TOM_AV2_XMIPP_CHECKSTACK checks a series of particles for use in xmipp-reconstructions 
%
%   stat=tom_av2_xmipp_checkstack(in_sel,exp_size)
%
%PARAMETERS
%
%  INPUT
%   in_sel           sel-file with particles to be checked
%   exp_size         size stack should be
%   exp_pix_range    ([mean+-10std]) number of std values should be
%   norm             ('none') ... only mean0+1std implemented           
%   subset           (all) length of random subset
%   verose           (1) verbose flag  
%
% OUTPUT
%   stat            pixel value statistic
%   error_count     abs number of errors (should be 0)
%
%EXAMPLE
%    st=tom_av2_xmipp_checkstack('23_corrf_high_128.sel',[128 128],5,'mean0+1std');
%    figure; plot(st.mean(:));
%    tom_dev(st.mean(:))
%    tom_dev(st.std(:))
%    hist(st.std(:))
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 03/02/10
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

try
    sel=importdata(in_sel);
catch Me
    if (verbose)
        disp(['Cannot read ' in_sel]);
        disp(Me.message);
    end;
    error_count=1;
    stat=-1;
    return;
end;

if (isempty(sel))
    if (verbose)
        disp(['Empty ' in_sel]);
    end;
    error_count=1;
    stat=-1;
    return;
end;

if (nargin < 3)
    exp_pix_range=10;
end;

if (nargin < 4)
   norm='none';
end;

if (nargin < 5 || isempty(subset) || strcmp(subset,'all') ) 
    subset=1:length(sel.textdata);
else
    tmp=randperm(length(sel.textdata));
    if (subset > length(tmp))
       if (verbose)
            disp(['subset(' num2str(subset) ') > ' in_sel '(' num2str(length(sel.textdata))  ')' ]);
            disp(['reducing subset to: ' num2str(length(sel.textdata))]);
        end;
        subset= length(tmp);
    end;
    subset=tmp(1:subset);
end;

if (nargin < 6)
   verbose=1;
end



all_mean=ones(length(subset),1).*-1;
all_min=ones(length(subset),1).*-1;
all_max=ones(length(subset),1).*-1;
all_std=ones(length(subset),1).*-1;
all_var=ones(length(subset),1).*-1;
all_sz=ones(length(subset),2).*-1;

mean=-1;
max=-1;
min=-1;
std=-1;
variance=-1;
not_read={};
zz=0;
error_count=0;
in.Value=zeros(exp_size);

for i=1:length(subset)
    
    f_name=sel.textdata{i};
    
    try
        in=tom_spiderread(f_name);
        [mean, max, min, std, variance] = tom_dev(in.Value,'noinfo');
    catch Me
        zz=zz+1;
        if (verbose)
            disp(['read error: ' f_name]);
        end;
        not_read{zz}=f_name;
        disp(Me.message);
        error_count=error_count+1;
     end;
    
    if (sum(size(in.Value)==exp_size)~=2)
        if (verbose)
            disp(['size error: ' f_name]);
            disp(num2str(size(in.Value)));
        end;
        error_count=error_count+1;
    end;
    
    if (std==0)
        if (verbose)
            disp(['std error: ' f_name]);
        end;
        error_count=error_count+1;
    end;
    
    if (isnan(mean))
        if (verbose)
            disp(['nan error: ' f_name]);
        end;
        error_count=error_count+1;
    end;
    
    if ((min >  ((mean-exp_pix_range*std) ))  &&  (max <  ((exp_pix_range*std) + mean)) ) ==0
        if (verbose)
            disp(['data range error: ' f_name]);
        end;
        error_count=error_count+1;
    end;
    
    if (strcmp(norm,'mean0+1std'))
        if (abs(mean) > 0.05 &&  abs(std-1) > 0.05  )
            if (verbose)
                disp(['norm error: ' f_name]);
            end;
            error_count=error_count+1;
        end;
    end;
    
    all_mean(i)=mean;
    all_max(i)=max;
    all_min(i)=min;
    all_std(i)=std;
    all_var(i)=variance;
    all_sz(i,:)=size(in.Value);
    if (mod(i,10000)==0 && verbose)
        disp(num2str(i));
    end;
end;

stat.mean=all_mean;
stat.min=all_min;
stat.max=all_max;
stat.std=all_std;
stat.var=all_var;
stat.sz=all_sz;
stat.not_read=not_read;

if (verbose)
    flInfo='info';
else
    flInfo='noinfo';
end;

me=tom_dev(all_mean,flInfo);
if (verbose)
    disp(['all_means: (mean) ' num2str(me)]);
end;
[me,ma,mi,std]=tom_dev(all_std,flInfo);
if (verbose)
    disp(['all_std: (mean) ' num2str(me)]);
    disp(['all_std: (std) ' num2str(std)]);
end;
[me,ma,mi,std]=tom_dev(all_sz,flInfo);
if (verbose)
    disp(['all_size: (mean) ' num2str(me)]);    
    disp(['all_size: (std) ' num2str(std)]);
end;

