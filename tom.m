%TOM toolbox to do tomography
%
%See also
%
%Acquisition
%TOM_TILT_SERIES      Tool to acquire a tilt series (WDN & FB) 
%
%Analysis
%TOM_CCC              calculates normalized 3d cross correlation coefficient (FF)
%TOM_CM               calculates the center of mass (1D, 2D, 3D, AL)
%TOM_CORR             computes 3d cross correlation function (FF)
%TOM_CREATE_CTF       calculates 2D or 3D CTF (pure phase contrast)(2D, 3D, FF)
%TOM_CTF              calculates and plots CTF (pure phase contrast)(1D, FF)
%TOM_CTF_ZERO         calculates zeros of CTF (pure phase contrast, FF)
%TOM_DEV              calculates deviation for images (2D, 3D, AF, tested FF)
%TOM_HIST3D           calculates a histogram of 3D data(SN)
%TOM_ORCD             performs three-dim orientational search (FF)
%TOM_MOI              calculates the moments of inertia (3D, FF)
%TOM_PEAK             determines the maximum/minimum value (1D, 2D, 3D, AL & GS, tested FF)
%TOM_XRAYCORRECT      corrects X-Ray defects on EM-Images (2D, AL)
%
%Average
%TOM_AVERAGE2D        tool for 2D averaging of single particles (WDN &FB)
%
%Display
%TOM_CTFFIT           tool to determine the defocus (2D, WDN)
%TOM_DSPCUB           visualization of a 3D volume in a gallery (3D, AF)
%TOM_EMBROWSE         browser for EM files (2D,3D, WDN)  
%TOM_IMAGESC          similar to imagesc, transposes the image (EM style)(WDN)
%TOM_INTERVOL         visualization tool for 3D data loaded in memory (3D, WDN)                               
%TOM_ISOSURFACE       tool for isosurface representation (SN)
%TOM_PARTICLES        tool for interactive 3D particles picking (3D, WDN, tested SN)
%TOM_PICKER           tool for interactive 3D particles picking loaded in memory (FF)
%TOM_VOLXY            visualization tool for 3D volumes (SN)
%
%Filtering & Transformation
%TOM_BANDPASS         performs bandpass filtering of image or volume (FF)
%TOM_CUT              filtering in Fourier Space (AL)
%TOM_FILTER           convolutes with spherical or quadratic kernel (FF)
%TOM_FOURRIER         calculates the fourier transformation of an image (1D, 2D, 3D, AL, tested FF)
%TOM_IFOURRIER        calculates the inverse fourier transformation of an image (1D, 2D, 3D, AL, tested FF)
%TOM_LAPLACE          performs Laplace filtering of image or volume (FF)
%TOM_OSCAR            performs 3D cross correlation of template and volume (FF)
%TOM_PS               calculates power spectrum of an image (2D, 3D, FF, tested FF)
%TOM_SHIFT            shifts image by a vector (even subpixels) (1D, 2D, 3D, FF)
%TOM_SMOOTH           performs smoothing of borders (2D, 3D, FF, tested FF)
%TOM_SMOOTHZ          performs smoothing of borders along z only (2D, 3D, FF, tested FF)
%TOM_WEDGE            calculates wedge - to be used as a filter for tomo (3D, FF)
%TOM_WEIGHT3D         weighting function for tomographic reconstruction (FF & FB)
%
%Geometrical shapes
%TOM_CIRCLE           generates a circle (AL)
%TOM_CYLINDER         creates an 3D array with a cylinder inside (SN)
%TOM_CYLINDERMASK     masks volume with (fuzzy) cylinder (3D, FF)
%TOM_ELLIPSEMASK      masks volume with (fuzzy) ellipse (3D, FF)
%TOM_SPHERE           creates volume with a sphere (FF)
%TOM_SPHEREMASK       masks volume with (fuzzy) sphere  (2D, 3D, FF)
%
%Input/Output functions
%TOM_CONVERT          converts images to a different file format (AL & GS)
%TOM_EMHEADER         Adds a EM-header structure to a matrix (SN)
%TOM_EMREAD           reads data in EM-file format (SN)
%TOM_EMREADC          reads data in EM-file format Developer version (SN)
%TOM_ENWRITE          writes data in EM-file format (SN)
%TOM_EMWRITEC         writes data in EM-file format Developer version (SN)
%TOM_ISEMFILE         checks if data is in EM-file format (SN)
%TOM_ISMRCFILE        checks if data is in MRC-file format (SN)
%TOM_MRCREADCLASSIC   reads data in MRC-file format (SN)
%TOM_MRCSTACK2EMSERIES  converts a MRC-file format stack in a series of single EM-files (SN)
%TOM_MRCSTACK2EMSTACK   converts a MRC-file format stack in a EM-format stack (SN)
%TOM_MRCWRITECLASSIC  writes data in MRC-file format (SN)
%TOM_PDB2EM           Converts PDB-files to EM files (SN)
%TOM_PDBREAD          reads a Protein Data Bank (pdb) file into a MATLAB structure (SN)
%TOM_RAWREAD          user interface for reading raw data (SN)
%TOM_READEMHEADER     reads only the header of an EM-file (SN)
%TOM_READMRCHEADER    reads only the header of an MRC-file (SN)
%
%Miscellaneous
%TOM_ERROR            generates Gaussian or Poisson noise (AL)
%TOM_MAKEMOVIE        generates an avi-file format movie file from a 3D stack (SN)
%
%Reconstruction
%TOM_BACKPROJ3D       performs 3D backprojection (SN)
%TOM_REC3D            a tool for tomographic 3D reconstructions (WDN)
%TOM_RECONSTRUCTION3D performs an aligned, weighted backprojection (SN & FB)
%TOM_RECPARTICLES     to reconstruct particles with high resolution, from tom_particles (3D, WD, tested SN)
%TOM_SETMARK          tool for the alignment of projections (2D, WDN, tested FF)
%
%Spatial Transformation
%TOM_BIN              performs binning of 1D, 2D or 3D images (1D,2D,3D , FF, tested SN)
%TOM_CART2CYL         transforms 3D-volumes from cartesian to cylinder coordinates (FF) 
%TOM_CART2POLAR       transforms 2D-images from cartesian to polar coordinates (2D, FF)
%TOM_CART2SPH         transforms 3D-volumes from cartesian to spherical coordinates (3D, FF)
%TOM_CYL2CART         transforms 3D volume from polar to cartesian coordinates(FF)
%TOM_LIMIT            sets limits to the values of an image (2D, 3D, AL)
%TOM_MIRROR           creates the mirror of an image (2D, 3D, AL)
%TOM_MOVE             moves the image in x,y,z direction (2D, 3D, AL, tested FF)
%TOM_PASTE            pastes array a into array b (AL)
%TOM_POINTROTATE      rotates a 3D vector point (FF)
%TOM_POLAR2CART       transforms 2D-images from polar to cartesian coordinates (2D, FF)
%TOM_RED              extracts a subimage/subvolume from the input image/volume (AL)
%TOM_ROTATE           rotates a 2D or 3D data  (3D, SN & FB)
%TOM_SPH2CART         transforms 3D-volumes from polar to cartesian coordinates (3D, FF)
%TOM_SYMREF           does n-fold symmetrization of a 3D reference (3D, FF)
%
%Utilities
%TOM_CHANGEMARKERFILE    tool to modify a marker file (WDN)
%TOM_CHOOSER             subjective classification of particles based on a motive list (3D, FF, tested FF)
%TOM_SORTSERIES          copy, correct x-rays and sort a tiltseries by angles (SN)
%TOM_TILT_LINES          plots marker points (FF) 
%
%To Do:
%tom_proj3d2d  
%tom_phantom3d
%tom_proj3d2d
%tom_mrc2em
%tom_mrcread             reads data in MRC-file format, FEI style (SN)
%tom_spiderread          reads data in SPIDER-file format (SN)
%tom_infoseries
%tom_imadj
%tom_symref2d            does n-fold symmetrization of a 2D reference (2D, FF)
%tom_tiltstrech
%tom_unbend
%tom_crowther            tool for estimating the resolution with the s.c. Crowther criteria (SN)
%tom_average3d
%tom_compare             compares 2 arrays in Fourier space (2D, 3D, MR)
%
%Installation
%makefun                 compiles Functions (please extend after creation of Functions ...)	
% 
%    Electron Tomography toolbox of the
%    Max-Planck-Institute for Biochemistry
%    Dept. Molecular Structural Biology
%    82152 Martinsried, Germany
%    http://www.biochem.mpg.de
% 
%    Copyright, 27/11/02
%    Last changed: 30/04/04
%
%
%
