function vol=tom_simulate_backproj_2d(log,st,l,b);
%TOM_SIMULATE_BACKPROJ_2D creates ...
%
%   vol=tom_simulate_backproj_2d(log,st,l,b);
%
%PARAMETERS
%
%  INPUT
%   log                 ...
%   st                  ...
%   l                   ...
%   b                   ...
%  
%  OUTPUT
%   vol                 ...
%
%EXAMPLE
%   ... = tom_simulate_backproj_2d(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

vol=zeros(size(st,1),size(st,2),size(st,1),'single');
st=single(st);

figure;
for i=1:size(log,1)
    tiltaxis=log(i,1);
    tiltangle=log(i,2);
    eu_all(i,:)=tom_sum_rotation([0 0 tiltangle; 270 90 tiltaxis],[0 0 0 ; 0 0 0]);
end;


for i=1:size(st,3);
    tiltaxis=log(i,1);
    tiltangle=log(i,2);
    proj=st(:,:,i);

    eu=tom_sum_rotation([0 0 tiltangle; 270 90 tiltaxis],[0 0 0 ; 0 0 0]);
    
    %generic sample thickness
    if (l==b)
        thickn=l;
    else
        thickn=l.*abs(cos((tiltaxis.*(pi./180)))) + b.*abs(sin((tiltaxis.*(pi./180)))); 
    end;
    
    [a,b,iso]=tom_calc_pixel_thresh(-tom_norm(proj,1),5000,0.01);
%     
    mask=zeros(size(proj));
    mask=tom_paste(mask,ones(120,50),[20 55]);
    mask=tom_filter(mask,6);
    %proj=mask.*-(iso~=0);
    proj=mask.*-(iso);
    
    proj=tom_smooth(proj,round(size(proj,1)./10));
    %proj=tom_weight3d_euler(proj,size(vol,1)./2-2,eu(1),eu(2),log(:,5),i,thickn);
    
     proj=tom_rotate(proj,90,'linear');
    w=tom_weight2d_exact([160 160],eu,eu_all,thickn);
    proj = real(ifft2((fft2(proj).*fftshift(w) ))); 
     proj=tom_rotate(proj,-90,'linear');
    
    tom_imagesc(double(proj)); drawnow;
    tom_backproj3d_euler(vol,proj,eu(1),eu(2),eu(3),[0 0 0]);
end;

disp('end');










