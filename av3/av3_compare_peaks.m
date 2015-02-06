function [peaklist, comp_motl] = av3_compare_peaks(motivelist, corrvolume, r)
%AV3_COMPARE_PEAKS compares peaks
%
%   [peaklist, comp_motl] = av3_compare_peaks(motivelist, corrvolume, r)
%
%   AV3_COMPARE_PEAKS reads out values of CORRVOLUME at given coordinates which
%   are defined in MOTIVELIST. In a sphere of radius R the local maximum is 
%   determined. The procedures can also read out the angles from a stack containing
%   the indexes of orientations given by ANGLES and add it to a new motivelist
%   COMP_MOTL.
%   ANGLES not implemented yet...
%Input
%   motivelist     : Motivelist of dim 20 x no_of_particles  (for more type help av3_collectparticles)
%   corrvolume     : Volume containing the correlation function of a particle to be compared with
%   r              : radius
%
%Output
%   peaklist       : matrix of dim 3 x no_of_particles. Column 1: part no; 2: CCC_in; 3: CCC2_out
%   comp_motl      : motivelist of compared particle
%
%   SEE ALSO
%   AV3_COLLECTPARTICLES, oscar
%
%   started 
%   08/10/02 FF
peaklist = zeros(3,size(motivelist,2));
comp_motl = 0*motivelist;
for i = 1:size(motivelist,2),
    x = motivelist(8,i); y = motivelist(9,i);z = motivelist(10,i);
    CCC1 = motivelist(1,i);
    box = tom_spheremask(corrvolume(x-r:x+r,y-r:y+r,z-r:z+r),r);  
    %box = box.*tom_sphere(r); % noch red(circ)
    [rs, CCC2] = tom_peak(box);
    x = rs(1) + x + r; y = rs(2) + y + r; z = rs(3) + z + r;
    comp_motl(8,i) = x;comp_motl(9,i) = y;comp_motl(10,i) = z;
    comp_motl(1,i) = CCC2;
    peaklist(1,i) = i;
    peaklist(2,i) = CCC1;
    peaklist(3,i) = CCC2;
end