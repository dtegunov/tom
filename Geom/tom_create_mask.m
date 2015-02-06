function mask_out=tom_create_mask(mask_st)
%TOM_CREATE_MASK creates a mask according to the given structure
%
%   mask_out=tom_create_mask(mask_st)
%
%PARAMETERS
%
%  INPUT
%   mask_st             mask structure:
%                           mask.Apply:  1 built mask according to the struct values, 2 use default values, 0 create ones
%                           mask.Value:  vector with parameters [size,radius,smooth,center ...]
%                           mask.Method: 'sphere'
%                           'sphere3d' 'cylinder3d' 'rectangle'  
%  
%  OUTPUT
%   mask_out            mask
%
%this is used for compatibility with tom_filtergui (AK)
%
%EXAMPLE
% mask_st.Apply=1; mask_st.Value=[256 256 30 10]; mask_st.Method='sphere';
% 
% mask=tom_create_mask(mask_st);
%
% figure; tom_imagesc(mask); 
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 25/01/06
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

tag=mask_st.Type(size(mask_st.Type,2)-1:size(mask_st.Type,2));

if (isempty(mask_st.Value.size_x) | isempty(mask_st.Value.size_y))
    error('size must be specified !!');
end



if (mask_st.Apply==0)
    if (strcmp(tag,'3d'))
        mask_out=ones(mask_st.Value.size_x,mask_st.Value.size_y,mask_st.size_z);
    else
        mask_out=ones(mask_st.Value.size_x,mask_st.Value.size_y);
    end;
    return;
end;



if (strcmp(tag,'3d'))
    sz=[mask_st.Value.size_x mask_st.Value.size_y mask_st.Value.size_z];
else
    sz=[mask_st.Value.size_x mask_st.Value.size_y];
end;



switch lower(mask_st.Type)

    case 'sphere2d'
        if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.radius=round(min(sz)./2);
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
        end;
        %check for empty Values
        if (isempty(mask_st.Value.radius))
            st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            st.Value.sigma=0;
        end;
        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
        end;
        mask_out=tom_spheremask(ones(sz),mask_st.Value.radius,mask_st.Value.sigma,[mask_st.Value.center_x mask_st.Value.center_y 1]);

      case 'sphere'
        if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.radius=round(min(sz)./2);
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
        end;
        %check for empty Values
        if (isempty(mask_st.Value.radius))
            st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            st.Value.sigma=0;
        end;
        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
        end;
        mask_out=tom_spheremask(ones(sz),mask_st.Value.radius,mask_st.Value.sigma,[mask_st.Value.center_x mask_st.Value.center_y 1]);

          
        
        
        
        
    case 'sphere3d'
        if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.radius=round(min(sz)./2);
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=round(min(sz)./2)+1;
        end;
        %check for empty Values
        if (isempty(mask_st.Value.radius))
            st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            st.Value.sigma=0;
        end;

        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y) |  isempty(mask_st.Value.center_z))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=round(min(sz)./2)+1;
        end;
        mask_out=tom_spheremask(ones(mask_st.Value.size_x,mask_st.Value.size_y,mask_st.Value.size_z),mask_st.Value.radius,mask_st.Value.sigma,[mask_st.Value.center_x mask_st.Value.center_y mask_st.Value.center_z]);



    case 'cylinder3d'
        if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.radius=round(min(sz)./2);
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=round(min(sz)./2)+1;
        end;
        %check for empty Values
        if (isempty(mask_st.Value.radius))
            st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            st.Value.sigma=0;
        end;
        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y) |  isempty(mask_st.Value.center_z))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=round(min(sz)./2)+1;
        end;
        mask_out=tom_cylindermask(ones(mask_st.Value.size_x,mask_st.Value.size_y,mask_st.Value.size_z),mask_st.Value.radius,mask_st.Value.sigma,[mask_st.Value.center_x mask_st.Value.center_y mask_st.Value.center_z]);

        
    case 'rectangle'
        
        if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.radius=round(min(sz)./2);
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
        end;
        %check for empty Values
        if (isempty(mask_st.Value.radius))
            st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            st.Value.sigma=0;
        end;
        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y) |  isempty(mask_st.Value.center_z))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
        end;
        mask_out=tom_spheremask(ones(sz),st.Value.radius,st.Value.sigma,st.Value.center);
        

    case 'ellipse'
        if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.mask_st.Value.radius1=round(min(sz)./2);
            mask_st.Value.mask_st.Value.radius2=round(min(sz)./2)-0.5.*round(min(sz)./2);
            mask_st.Value.mask_st.Value.radius3=1;
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=1;
        end;
        %check for empty Values
        if (isempty(mask_st.Value.radius1))
            mask_st.Value.radius1=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.radius2))
            mask_st.Value.radius2=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.radius3))
            mask_st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            mask_st.Value.sigma=0;
        end;
        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y) |  isempty(mask_st.Value.center_z))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=1;
        end;
        mask_st.Value.radius3=1;
        mask_out=tom_ellipsemask(ones(sz),mask_st.Value.radius1,mask_st.Value.radius2,mask_st.Value.radius3,mask_st.Value.sigma,[mask_st.Value.center_x mask_st.Value.center_y mask_st.Value.center_z]);
        
        
    case 'ellipse3d'
         if (mask_st.Apply==2)
            %Default Settings
            mask_st.Value.mask_st.Value.radius1=round(min(sz)./2);
            mask_st.Value.mask_st.Value.radius2=round(min(sz)./2)-0.5.*round(min(sz)./2);
            mask_st.Value.mask_st.Value.radius3=round(min(sz)./2);
            mask_st.Value.sigma=round(min(sz)./20);
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=round(min(sz)./2)+1;
        end;
         %check for empty Values
        if (isempty(mask_st.Value.radius1))
            mask_st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.radius2))
            mask_st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.radius3))
            mask_st.Value.radius=round(min(sz)./2);
        end;
        if (isempty(mask_st.Value.sigma))
            mask_st.Value.sigma=0;
        end;
        if (isempty(mask_st.Value.center_x) | isempty(mask_st.Value.center_y) |  isempty(mask_st.Value.center_z))
            mask_st.Value.center_x=round(min(sz)./2)+1;
            mask_st.Value.center_y=round(min(sz)./2)+1;
            mask_st.Value.center_z=round(min(sz)./2)+1;
        end;
        mask_out=tom_tom_ellipsemask(ones(sz),mask_st.Value.radius1,mask_st.Value.radius2,mask_st.Value.radius3,mask_st.Value.sigma,mask_st.Value.center);

    case 'roipoly'
        mask_out = ones(mask_st.Value.size_x,mask_st.Value.size_y);
        pm = zeros(mask_st.Value.size_x,mask_st.Value.size_y);
        for k=1:size(mask_st.Polygons,2)
            coord = mask_st.Polygons{k}./2^(mask_st.Value.binning);
            pm = pm + roipoly(mask_out,coord(:,2)',coord(:,1)');
        end
        pm(find(pm>1))=1;
        mask_out = -pm+1;
        
end;

%check for rotation of mask
if (isfield(mask_st.Value,'angle_phi') && isfield(mask_st.Value,'angle_psi') && isfield(mask_st.Value,'angle_theta'))

    if (not(isempty(mask_st.Value.angle_phi)) || not(isempty(mask_st.Value.angle_phi)) || not(isempty(mask_st.Value.angle_theta)))

        if (isempty(mask_st.Value.angle_phi))
            mask_st.Value.angle_phi=0;
        end;

        if (isempty(mask_st.Value.angle_psi))
            mask_st.Value.angle_psi=0;
        end;

        if(isempty(mask_st.Value.angle_theta))
            mask_st.Value.angle_theta=0;
        end;

        angle=[mask_st.Value.angle_phi mask_st.Value.angle_psi mask_st.Value.angle_theta];
        mask_out=tom_rotate(mask_out,angle);
    end;
end;




