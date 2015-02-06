% ALIGNMENT
%   AV3_TRANS_ROT_ALIG   aligns 3D subtomograms to reference and averages
%                           iteratively
%   AV3_SCANANGLES_EXACT aligns 3D subtomograms to reference and averages
%                           iteratively (same as AV3_TRANS-ROT-ALIG)
%   AV3_TRANS_ROT_ALIG_NOPHI aligns 3D subtomograms to reference
%   AV3_GLOBAL_SCAN      aligns 3D subtomograms to reference - exhaustive 
%                           scan!
%   AV3_PHIALIG          aligns 3D subtomograms along PHI
%   AV3_AVERAGE_EXACT    averages subtomograms
%   AV3_AVERAGE          averages 3D-particles from entire TOMOGRAM
%   AV3_AVERAGE_WEI      performs WEIGTED averaging of 3D-particles from 
%                           entire tomograms 
%   AV3_CREATEMOTL       creates motivelist according to MOLMATCH run
%   AV3_DOCK             aligns 3D reference (e.g. from X-ray) to density
%                           (from EM)
%   AV3_DOCK_NONORM      just as AV3_DOCK but different normalization
%   AV3_DOCK_PHASE_ONLY  just as AV3_DOCK but different normalization
%   AV3_INDEX2ANGLE      creates Euler angles from molmatch indices
%   AV3_MOTLANALYZE      plots convergence of MOTL
%
% REC_DATA
%   AV3_EXTRACT_PROJS    projects subtomograms
%   AV3_RECPARTICLES     reconstructs subtomograms according to MOTL
%   AV3_FASTRECPARTICLES reconstructs subtomograms according to MOTL
%   AV3_LOCATEONPROJ     computes 3D coordinates to 2D coordinates in
%                           projection
%
% UTILS
%   AV3_BANDPASS         performs bandpass filtering of image or volume
%   AV3_CART2SPH         transforms Cartesian to spherical coordinates
%   AV3_CYLSIM           can be used to symmetrize a volume cylindrically
%   AV3_LAGUERRE         creates Laguerre Polynomials
%   AV3_MATRIX2ANGLE     converts rotation matrix to angles
%   AV3_ROTAVERAGE       rotationally averages 3D volume
%   AV3_VOL2GALLERY      converts a 3D volume into a 2D gallery 
%   AV3_WEDGE            produces a wedge shaped array
%   ELLIPSEDET           determines ellipse parameters from parameters
%   ELLIPSEPOS2ANGLE     determines Euler angles from ellipse parameters
%   NORMALVEC            determines Euler angles of particles on a sphere
%   GAU                  error of data and Gaussian (for fitting)
%   GAUGAU               error of data and Double Gaussian
%   PEAK_DET_2           Determination of the maximum value
%   TOM_ZOOM             zooms
%
% VISU
%   AV3_CHECK...         displays aligned particles
%   ORIENTATIONS 
%   AV3_MOTL2MASK        creates a mask from a MOTL
%   AV3_VISUALIZE        visualizes 3D-particles (similar to AV3_CHECK)
%   TOM_CHOOSER          interactively "classifys" particles in tomogram
%
% ANALYSIS AND CLASS
%   AV3_RESDET           determines resolution of average
%   AV3_VARIANCE         computes variance of average
%   AV3_SYMTEST          performs rotational symmetry test
%   AV3_COMPNORMALS      computes particles' normalvectors to determined
%                           orientations
%   AV3_CONSTR_CCC       computes correlation of particles and template
%   CTR_CLASS determines centers of masses
%   TOM_RADPS            computes radially averaged powerspectrum
%   AV3_EIGFAC2VOL       converts eigenfactors to eigenvectors
%   AV3_MOIS             computes principal moments of inertia of particles
%   AV3_CCCMAT           computes correlation- and covariance matrix of particles
%   AV3_PREP_CLASS_IN_EM perpares stacks for classification in EM
%   
%
% ------------------------------------------------------------------------
%   Reference for this toolbox:
%   Förster, F., Medalia, O., Zauberman, N., Baumeister, W. & Fass, D.
%   "Retrovirus Envelope Protein Complex Structure in Situ Determined by 
%   Cryo-Electron Tomography". Proc Natl Acad Sci U S A 102, 4729-4734 (2005).
% ------------------------------------------------------------------------
%
%   Copyright (c) 2005
%   Max-Planck-Institute for Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom
%
% last update 03/31/05 FF


%   AV3_AVERAGE          performs averaging of 3D-particles - particles in whole tomogram
%   AV3_AVERAGE2         performs averaging of 3D-particles - particles stored in subtomograms
%   AV3_CART2SPH         Transform Cartesian to spherical coordinates - single points!
%   AV3_COLLECTPARTICLES collects particles from ccf, angles and volume
%   AV3_COMPARE_PEAKS    compares peaks - use to compare one X-CCF with another
%   AV3_CREATEMOTL       creates MOTL from CCF
%   AV3_CUTVOLUME        cuts large volume into pieces with overlap
%   AV3_FILLMOTL         creates MOTL from PEAKLIST
%   AV3_FINESCAN	       fragment ... delete?
%   AV3_INDEX2ANGLE      creates Euler angles from oscar indices
%   AV3_LAGUERRE         creates Laguerre Polynomials (for H-atom)
%   AV3_MERGEVOLUME      merges subvolume to large (CCF) volume
%   AV3_MOTLANALYZE      plots convergence of MOTL
%   AV3_PARTICLESELECT   M-file for av3_particleselect --- ?? tom_picker
%   AV3_PROJANGLES       calculates projection angles of individual projections -
%   NORMALVEC	           determine normal of points on a sphere
%   AV3_SCANANGLES_MOVINGMASK same as AV3_SCANANGLES_EXACT but with local
%                       normalization
% AV3_DOCK aligns 3D reference (e.g. from X-ray) to particle (from EM)
% ELLIPSEDET           determine ellipse parameters from parameters
% CREATE_REF           example for creating EM-reference from pdb-file
% PEAK_DET2            TOM_PEAK including interpolation - from Gabi
%			to be finished
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
