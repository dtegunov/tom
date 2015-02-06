function out=tom_rotate(varargin)
%TOM_ROTATE performs a 2d or 3d rotation, depending on the input
%
%   out=tom_rotate2(varargin)
%   out=tom_rotate(in,[euler_angles],interp,[center])
%
%NOTE
%   This function works as an interface. The actual computation is done
%   in the C-Function tom_rotatec
%
%PARAMETERS
%
%  INPUT
%   in                  image or Volume as single!
%   euler_angles        rotation matrixs [phi, psi, the] (rot X1 X3 Z2)
%   interp              interpolation only 'linear' implemented
%   center(optional)    in=image  [centerX centerY]
%                       in=Volume [centerX centerY centerZ]
%                       
%                       WARNING: Convention for rotation center
%                       Matlab / TOM : Because Matlab starts counting array elements at 1,
%                       rotation center will be size / 2 +1 by default
%
%                       C++ / Python: That equals size / 2 for real
%                       computer languages which start counting elements at 0 
%
%                       ( => Beware, Matlab is no real computer language. Is this a feature? )
%
%
%   taper               flag for tapering before roation
%  
%  OUTPUT
%   out                 Image rotated
%
%ROTATION INFORMATION
% 
% tom_rotatec.c uses the following 3d Euler rotation 
% (first index: line, second index: column):
% 
% rm00=cospsi*cosphi-costheta*sinpsi*sinphi;
% rm10=sinpsi*cosphi+costheta*cospsi*sinphi;
% rm20=sintheta*sinphi;
% rm01=-cospsi*sinphi-costheta*sinpsi*cosphi;
% rm11=-sinpsi*sinphi+costheta*cospsi*cosphi;
% rm21=sintheta*cosphi;
% rm02=sintheta*sinpsi;
% rm12=-sintheta*cospsi;
% rm22=costheta;
% 
%  This is a 3d Euler rotation around (zxz)-axes
%  with Euler angles (psi,theta,phi)
%  
%  For more information about 3d Euler rotation visit:
%  wiki -> Definitions -> Overview 3d Euler rotation
%  
%  To produce the same result as in EM use nx/2, ny/2
%  and nz/2, not !!!! +1
%
%       X-axis  corresponds to  phi=0     psi=0   the=alpha
%       Y-axis  corresponds to  phi=270   psi=90  the=alpha
%       Z-axis  corresponds to  phi=alpha psi=0   the=0
%
%EXAMPLE
%   im=tom_emread('pyrodictium_1.em');
%   out=tom_rotate(single(im.Value),40,'linear',[256 256],'taper');
%   
%   im=tom_emread('testvol.em');
%   im=single(im.Value);
%   out=tom_rotate(im,[30 20 95],'linear');
%
%REFERENCES
%
%SEE ALSO
%   TOM_VOLXYZ, TOM_REC3D, TOM_RECPARTICLES
%
%   created by FB 12/04/04
%   updated by FB 07/10/05
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

switch nargin
    case 5
        taper= varargin{5};
        center = varargin{4};
        ip=varargin{3};
    case 4,
        if (isnumeric(varargin{4})==0)
            taper=varargin{4};
            center = [ceil(size(varargin{1})./2)];
        else
            taper='no';
            center = varargin{4};
        end;
        ip=varargin{3};
    case 3,
        ip=varargin{3};
        center = [ceil(size(varargin{1})./2)];%bug fixed for uneven dims - FF
        taper='no';
    case 2,
        center = [ceil(size(varargin{1})./2)];%bug fixed FF
        ip = 'linear';
        taper='no';
    otherwise
        disp('wrong number of Arguments');
        out=-1; return;
end;
%parse inputs
in = varargin{1};
euler_angles=varargin{2};

if (strcmp(taper,'taper')==0);


    % allocate some momory
    out = single(zeros(size(in)));

    % call C-Function to do the calculations
    tom_rotatec(single(in),out,single([euler_angles]),ip,single([center]));

    out = double(out);

else

    cent=round((size(in)./2))+1;

    if (find(((center==cent)==0),1) )
        error(['taper only works with center n/2+1 center must be ' num2str(cent) ]);
    end;
        
        
    old_sz=size(in);
    sz=sort(size(in));
    max_sz=round(sz(size(sz,2)).*sqrt(2));

    if (size(sz,2)==2)
        max_size=[max_sz max_sz];
    else
        max_size=[max_sz max_sz max_sz];
    end;

    
    in=tom_taper(in,max_size);
    
    % allocate some momory
    out = single(zeros(size(in)));

    % calc new center
    center = [ceil(size(in)./2)];
    
    
    % call C-Function to do the calculations
    tom_rotatec(single(in),out,single([euler_angles]),ip,single([center]));
    
    out=tom_cut_out(out,(round((max_size-old_sz)./2)+1),old_sz,'nofill');
    out = double(out);

end;

