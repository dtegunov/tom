function stack = ctr_class(motl,particlefilename,threshold,filter)
% CTR_CLASS determines centers of masses
%
%   stack = ctr_class(motl,particlefilename,threshold,filter)
%
%   Procedure was developed for classifion of central transporter in NPC
%   project. A classification according to the centers of mass was
%   performed in this case. CTR_CLASS determines the centers of mass for
%   each NPC-CTR. 
%   The procedure inverts the gray values prior to computing the center of
%   mass (dark ~ dense), the THRESHOLD is chosen as the origin of mass
%   scale. 
%   Attention: In this version there are still hardcoded parameters.
%
% PARAMETERS
%  INPUT
%   motl                motivelist
%   particlefilename    filename of 3D-particles to be aligned and averaged
%   threshold           origin of mass scale - values above this value are
%                           set to zero. should be the grey value of the
%                           desired reference density (e.g. solvent)
%   filt                parameter for the filtering of particles - if 0 no
%                           filter    
% 
%  OUTPUT
%   stack               1st column: contains absolute mass
%                       2nd-4th: [x y z] coordinates of center of mass 
%                       5th: occupied voxels
%
% FF
% modified MB 23/10/03
%   last change 04/01/05 FF - updated docu
mask = tom_ellipsemask(ones(64,64,64),9,9,11,0,[33 33 37]);
npix = sum(sum(sum(mask)));
for ind = 1:size(motl,2)
    if ( (motl(20,ind) == 1) )
        disp(['now particle # ' num2str(ind) ' ...']);
        psi_opt = motl(18,ind);
        phi_opt = motl(17,ind);
        the_opt = motl(19,ind);
        tshift(1) = motl(11,ind);
        tshift(2) = motl(12,ind);
        tshift(3) = motl(13,ind);
        name = [particlefilename '_' num2str(ind) '.em'];
        particle = tom_emread(name);
        if (filter > 0) 
            particle.Value = tom_filter(particle.Value,filter);
        end;
        particle = double(tom_rotate3d(tom_shift(particle.Value,-tshift),-psi_opt,-phi_opt,-the_opt));
        particle = tom_limit(-particle+threshold,0,3,'z');%thresholden und invertieren
        particle = mask.*particle;
        figure(1);tom_dspcub(particle);
        mass = sum(sum(sum(particle)));
        r = tom_cm(particle);
        binparticle = zeros(64,64,64);
        indx = find(particle > 0.0000001);
        binparticle(indx) = 1;
        nels = sum(sum(sum(binparticle)));
        figure(2);tom_dspcub(binparticle);
        stack(1,ind)=mass;
        stack(2,ind)=r(1);
        stack(3,ind)=r(2);
        stack(4,ind)=r(3);
        stack(5,ind)=nels;
    end;
end;