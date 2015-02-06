function particleStack = tom_os3_generateParticleStack(pickList,job,norm,img,binning)
%tom_os3_generateParticleStack
%   
% 
% 
%   tom_os3_generateParticleStack(pickList,job,norm)
%
%PARAMETERS
%
%  INPUT
%       pickList
%       job
%       norm
%       binning
%  OUTPUT
%       
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by TH2 07/07/07
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



%% load the img
if(nargin < 4)
    img = tom_emreadc3(job.volumeFile);
    img = img.Value;
end;
templateSize = job.templateSize;

%% allocate memory for the particles
%particleStack = zeros(templateSize(1),templateSize(2),length(pickList),'single');
particleStack = [];
mask = tom_os3_sphereMask(zeros(templateSize,'single'));


if(job.dimension == 2)
    %%  for each hit in the picklist
    for i =1:length(pickList)

        pick = pickList{i};
        x = pick.coordinates(1);
        y = pick.coordinates(2);
        %cut out the particle and align it roughly, norm it according to norm
        %flag
        part = tom_cut_out(img,[x-templateSize(1)/2 y-templateSize(2)/2],[templateSize(1) templateSize(2)],'no-fill');
    %     part = tom_rotate(part,-job.angleListOrig(pick.angle));
        if(strcmp(norm,'phase') || strcmp(norm,'2std') || strcmp(norm,'3std') || strcmp(norm,'mean0+1std'))
            part = tom_rotate(tom_norm(part,norm),-job.angleListOrig(pick.angle));
%        else
 %           part = tom_rotate(part,-job.angleListOrig(pick.angle));
        end;
        %mask the edges out
        %part = part.*mask;
        particleStack(:,:,i) = tom_bin(part,job.options.modifications.binning);
    end;
else
    %%  for each hit in the picklist
    for i =1:length(pickList)

        pick = pickList{i};
        x = pick.coordinates(1);
        y = pick.coordinates(2);
        z = pick.coordinates(3);
        %cut out the particle and align it roughly, norm it according to norm
        %flag
        part = tom_cut_out(img,[x-templateSize(1)/2 y-templateSize(2)/2 z-templateSize(3)/2],[templateSize(1) templateSize(2) templateSize(3)],'no-fill');

        angle = job.angleListOrig{pick.angle};
        phi = angle(1);
        angle(1) = angle(2);
        angle(2) = phi;
        angle = -1 * angle;
        if(strcmp(norm,'phase') || strcmp(norm,'2std') || strcmp(norm,'3std') || strcmp(norm,'mean0+1std'))
            

            part = tom_rotate(tom_norm(part,norm),angle);
        else
            part = tom_rotate(part,angle);
        end;
        %mask the edges out
        part = part.*mask;
        particleStack(:,:,:,i) = tom_bin(part,binning);
    end;
end;
