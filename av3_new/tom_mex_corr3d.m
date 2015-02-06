%TOM_MEX_CORR32 does normalised cross correlation of a list of templates 
% with a list of particles, for a given set of rotations.
%
%   [out] = tom_mex_corr3d([type], size, templates, particles, angles,
%                mask, wedge_ref, wedge_part, type, [...])
%
%PARAMETERS
%
%  INPUT
%   type                Optional parameter specifying the data-type 
%                       with which the computatio is done. It can be 
%                       either 'double' or 'single'.
%   size                The size of the volumes. All templates, references
%                       wedges etc. must have this size. It is a 3-vector.
%   templates           The templates to do the correlation. If can be
%                       either: - A volume of dimension size.
%                               - A filename of an em-file to a volume.
%                               - A cell-array of either volumes or
%                                 filenames (may be mixed).
%   particles           The list of particles. In the same form as
%                       templates.
%   angles              A volume containing the rotaion angles. It must
%                       be a 3xN matrix with the [phi,psi,theta]-euler
%                       angles in degrees.
%   mask                A mask applied to the particles and the templates.
%                       It can - [] for no mask
%                              - a number to make a sphere with a given
%                                radius (sigma==0), as in TOM_SPHEREMASK
%                              - A volume of dimension size. 
%                       This mask is applied to the particles and the 
%                       rotated templates.
%   wedge_ref           The wedge of the templates. It can be either
%                       - [] to specity that no template has a missing 
%                         wedge.
%                       - A number specifying the opening angle of the
%                         missing wedge as in TOM_WEDGE.
%                       - A real volume of dimension size which is the weighting
%                         function of the frequences in fourier space (i.e.
%                         the missing wedge).
%                       - A cell array of the above. Than each template has
%                         its own, individual missing wedge. Otherwise they
%                         all have the same one.
%   wedge_part          The wedge of the particles in the same form as
%                       wedge_ref
%   type                Specifies what the function should do with the
%                       correlation volume. Currently there are two options
%                       implemented: 
%     'peak': Finds the peak of the correlation volume.
%             This option allows a 9th parameter which is a
%             mask where a peak is only regarded if the
%             corresponding voxel of the mask is ~= 0. Similar
%             to the parameter mask_ccf of TOM_AV3_ALIGN.
%     'cc'    Returns a cell array with each correlation volume.
%             Can have an optional flag which decides to return the
%             fft-shifted correlation volume or the original one. Omitting
%             or setting to [] defaults to true (i.e. do the fftshift).
%     'cc2em': Each correlation volume is saved as an
%             em-file. This option needs a 9th parameter. It
%             is a NxP cell array of strings (N number of
%             templates, P number of particles) where for each
%             combination of template and particle the
%             directory name is specified where to save the
%             em-file. The file-name then is composed of the
%             current index of templates, particles and
%             angles and is returned by the function. An optional 10th
%             parameter is a boolean value which specifies whether to do a
%             fftshift of the correlation value before saving it to file.
%             Omitting this parameter or setting to [] defaults to true
%             (i.e. do fftshift of the ccvolume).
%  OUTPUT
%    depends on the type.
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_ALIGN
%
%   created by Thomas Haller 12/17/07
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

