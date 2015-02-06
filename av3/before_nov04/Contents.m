% BASICS
%
% AV3_AVERAGE          performs averaging of 3D-particles - particles in whole tomogram
% AV3_AVERAGE2         performs averaging of 3D-particles - particles stored in subtomograms
% AV3_CART2SPH         Transform Cartesian to spherical coordinates - single points!
% AV3_COLLECTPARTICLES collects particles from ccf, angles and volume
% AV3_COMPARE_PEAKS    compares peaks - use to compare one X-CCF with another
% AV3_CREATEMOTL       creates MOTL from CCF
% AV3_CUTVOLUME        cuts large volume into pieces with overlap
% AV3_FILLMOTL         creates MOTL from PEAKLIST
% AV3_FINESCAN	       fragment ... delete?
% AV3_INDEX2ANGLE      creates Euler angles from oscar indices
% AV3_LAGUERRE         creates Laguerre Polynomials (for H-atom)
% AV3_MERGEVOLUME      merges subvolume to large (CCF) volume
% AV3_MOTLANALYZE      plots convergence of MOTL
% AV3_PARTICLESELECT   M-file for av3_particleselect --- ?? tom_picker
% AV3_PROJANGLES       calculates projection angles of individual projections -
% NORMALVEC	           determine normal of points on a sphere
% AV3_SCANANGLES_EXACT aligns 3D subtomograms to reference and averages
%                       iteratively
% AV3_SCANANGLES_MOVINGMASK same as AV3_SCANANGLES_EXACT but with local
%                       normalization
% AV3_DOCK aligns 3D reference (e.g. from X-ray) to particle (from EM)
% ELLIPSEDET           determine ellipse parameters from parameters
% CREATE_REF           example for creating EM-reference from pdb-file
% PEAK_DET2            TOM_PEAK including interpolation - from Gabi
%			to be finished
% 
% 
%-------------------------------------------------------------------------------------
%
% SPHHARM - sphercoord
%
% sp3_ylmexpand      expands volume (in spherical coordinates) into spherical
%                    harmonics, i.e. expansion coefficients are determined
% sp3_ylmexpandfft   the same but using fft - much faster - recommended
% sp3_ylmrec         performs summation of spherical harmonics - 'reconstruction'
%                    or inverse transformation (in spherical coordinates)
% sp3_ylmrecfft      the same with correct normalization for FFT version - recommended
% sp3_ylmrec_cartes  reconstructs on cartesian grid - interpolation in r is performed
% sp3_rotate         performs rotation in l,m space
% sp3_setuprotate    subfunction of sp3_rotate - sets up matrices for rot
% sp3_rotate2        work in progress
% sp3_orcd           orcd using spherical harmonics
% sp3_mexchangemat   calculates exchange matrix of coefficients - needed for sp3_orcd
% sp3_creafilter     creates filter for template matching approach
% sp3_testfilter     test for sp3_creafilter
%--------------------------------------------------------------------------------------
% SPHHARM - 2D
%
% template matching using steerable filters
%--------------------------------------------------------------------------------------
%
% H-BASIS
%
%package for H-atom basis functions
%--------------------------------------------------------------------------------------
%
% last update 
% FF 08/18/03
%
