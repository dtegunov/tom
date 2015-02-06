function [profile minima]=tom_get_sym_profile(volume,nsym,profile,profile_param,demo)
%tom_get_sym_profile calculates symmetry profile 
%
%   [profile minima]=tom_get_sym_profile(volume,nsym,profile,profile_param,demo)
%
%  tom_get_sym_axis calculates symmetry profile of a volume by applying
%  symmetry and projecting
%  
%  
%  
%
%PARAMETERS
%
%  INPUT
%   volume             input volume
%   nsym               symmetry
%   profile            flag for extracting the profile
%                           2proj,cent_proj    
%   profile_param      parameter for cent_proj
%   demo               demo mode (0/1)   
%
%  OUTPUT
%   profile            profile 
%   maxima             vector of maxima 
%
%
%EXAMPLE
%    
% %project 2 times in x direction 2 get pixel profile
% [profile minima]=tom_get_sym_profile(vol.Value,6,'2proj',0,1);
% 
% %use cent slice and average over profile_param from the middle
% [profile minima]=tom_get_sym_profile(vol.Value,6,'cent_proj',10,1);
%
%NOTE:
% 
% symmetry axis has to be y-axis!! (zxz)
% 
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by fb 
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

if (nargin<5)
    demo=0;
end;

vol_mid=floor(size(volume)./2)+1;
sz_vol=size(volume);

sym=tom_symref(volume,nsym);

%create pixel profile
if (strcmp(profile,'cent_proj'))
    slice=squeeze(sym(vol_mid(1),:,:));
    fun=sum(-slice(vol_mid-profile_param:vol_mid+profile_param,:),1);
end;

if (strcmp(profile,'2proj'))
    slice=squeeze(sum(sym,1));
    fun=-squeeze(sum(squeeze(sum(sym,1)),1));
end;


%extract maxima
d_fun=diff(fun);
%d_fun=gradient(fun);
[index b]=tom_crossing(d_fun,[],0);
zz=1;
for i=2:size(index,2)
    if (sign(d_fun(index(i)-1)) ==1 &&  sign(d_fun(index(i)+1))==-1 )
        index_max(zz)=index(i);
        zz=zz+1;
    end;
end; 
index_max=index_max+1;


if (demo==1)
    figure; plot(fun,'k','LineWidth',2);
    hold on; plot(index_max,fun(index_max),'ro'); hold off;
    figure; tom_imagesc(tom_rotate(slice,-90),'noinfo'); axis off;
    for i=1:length(index_max)
        hold on;  plot(repmat(index_max(i),sz_vol(1),1),1:sz_vol(1)); hold off;
    end;
end;


minima=index_max;
%function was inverted !!





%take 3rd and 5th maximum to cut out ATPase

%old staff
% index=index+1;
% 
% [ind_start start_val]=tom_nearestpoint(vol_mid(1),index_max);  




% stop=index_max(ind_start-2);
% start=index_max(ind_start-4)+2;
% 
% start2=index_max(ind_start+2);
% stop2=index_max(ind_start+4);
% 
% f=fun;
% 
% 
% if (demo==1)
%     figure; plot(f,'k','LineWidth',2);
%     axis off;
%     hold on; plot(index_max,f(index_max),'ro'); hold off;
%     hold on; plot(start_val,f(start_val),'go'); hold off;
%     figure; tom_imagesc(tom_rotate(slice,-90),'noinfo'); axis off;
%     
%     hold on; plot(repmat(start,sz_vol(1),1),1:sz_vol(1)); plot(repmat(stop,sz_vol(1),1),1:sz_vol(1)); 
%     plot(repmat(start2,sz_vol(1),1),1:sz_vol(1)); plot(repmat(stop2,sz_vol(1),1),1:sz_vol(1));  hold off;
% end;



