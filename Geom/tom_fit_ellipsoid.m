function [e_mask,param]=tom_fit_ellipsoid(vol,v_x,v_y,v_z,thick)
%  tom_fit_ellipsoid fits a ellipsoid to a given volume
%  
%     [e_mask,param]=tom_fit_ellipsoid(vol,v_x,v_y,v_z,thick)
%  
%  PARAMETERS
%  
%    INPUT
%     vol                 search volume
%     v_x                 search range 4 x radius 
%     v_x                 search range 4 y radius 
%     v_z                 search range 4 z radius 
%     thick               wall thickness                
%            
%
%
%    OUTPUT
%     e_mask              fitted ellipsoid mask
%     param               fitted parameter  
%                         param(1)=rx
%                         param(2)=ry
%                         param(3)=rz
%                         param(4)=posx
%                         param(5)=posy
%                         param(6)=posz
%
%
%  EXAMPLE
%   
%  matlabpool close;
%  matlabpool open local 8; 
% 
%  dummy_ellipse=tom_move(tom_ellipsemask(ones(64,64,64),21,11,6,0,[30 28 31])-tom_ellipsemask(ones(64,64,64),19,9,4,0,[30 28 31]),[2 2 2]);
%  [e_mask,param]=tom_fit_ellipsoid(dummy_ellipse,[5:5:30],[5:5:30],[10:5:30],2);
%  figure; tom_dspcub(e_mask); set(gcf,'Name','Fitted Ellipse');
%  figure; tom_dspcub(dummy_ellipse); set(gcf,'Name','Input Ellipse');
%
%   
%  NOTE
%  
%   Rotation is not fitted !!
%   ==> Input Vol main axes have 2 be parallel 2 cosy
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_ellipsemask
%  
%     fb 03/11/2011
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom


zz=0;
for ix=v_x
    for iy=v_y
       for iz=v_z
          zz=zz+1;          
          tmp_list{zz}=[ix iy iz];              
       end;
    end;
end;

tmp_v=ones(size(vol));
mid=floor(size(vol)./2)+1; 

for i=1:length(tmp_list)
    p_big=tmp_list{i}+round(thick./2);
    p_small=tmp_list{i}-round(thick./2);
    mask_big=tom_ellipsemask(tmp_v,p_big(1),p_big(2),p_big(3));
    mask_small=tom_ellipsemask(tmp_v,p_small(1),p_small(2),p_small(3));
    mask_fin=mask_big-mask_small;
    if (var(mask_fin(:))==0)
        pos(:,i)=[0 0 0];
        val(i)=-1;
        continue;
    end;
    if (i==212)
        disp(' ');
    end;
    
    cc=tom_corr(mask_fin,vol,'norm');
    [tmp_pos tmp_val]=tom_peak(cc);
    pos(:,i)=tmp_pos(1,:); 
    val(i)=tmp_val(1);
end;

[m_val m_pos]=max(val);

p_big=tmp_list{m_pos}+round(thick./2);
p_small=tmp_list{m_pos}-round(thick./2);
mask_big=tom_ellipsemask(tmp_v,p_big(1),p_big(2),p_big(3));
mask_small=tom_ellipsemask(tmp_v,p_small(1),p_small(2),p_small(3));
e_mask=mask_big-mask_small;
e_mask=tom_shift(e_mask,(pos(:,m_pos)-mid').*-1);
param=tmp_list{m_pos};
param(4)=pos(1,m_pos);
param(5)=pos(2,m_pos);
param(6)=pos(3,m_pos);

