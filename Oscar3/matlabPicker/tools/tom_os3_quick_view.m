function tom_os3_quick_view(align2d,filter,image,filename)
%tom_os3_quickview paste one array in another one
%
%   tom_os3_quickview(picklist,image)
%
%
%  
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   a                   ...
%
%EXAMPLE
%
% use for db
% tom_os3_quick_view('align2d_db.mat','','','/fs/sun17/lv01/pool/pool-nickell/tmp/scratch_fb/test_picking/pick/neagative_stain_20S_output/volumes2/20S_11.em');
%
%REFERENCES
%
%SEE ALSO
%   TOM_MOVE, TOM_PEAK, TOM_RED
%
%   created by AL 08/15/02
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

if nargin > 3
     if (ischar(align2d))
        load(align2d);
    end;
    align2d=align2d.get_subset(filename);
    image=tom_emread(align2d(1,1).filename);
    image=image.Value;
end;

if nargin < 3
    if (ischar(align2d))
        load(align2d);
    end;
    image=tom_emreadc3(align2d(1,1).filename);
    image=image.Value;
end;

if nargin < 2
    filter='NoFilter';
end;

if (ischar(image) && isempty(image)==0)
    image=tom_emread(image);
end;  

if (ischar(filter)==0 && isempty(filter)==0)
    image=tom_filter(image,filter);
end;

image=single(image);

figure; tom_imagesc(image); 


for i=1:size(align2d,2)
    hold on; plot(align2d(1,i).position.x,align2d(1,i).position.y,'ro'); hold off;
end;
    
[t fname ext]=fileparts(align2d(1,1).filename);
set(gcf,'Name',[num2str(i) ' Particles on ' fname ext ]);

disp('');



