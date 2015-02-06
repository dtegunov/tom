% TOM Toolbox
% Version 1.99
% 
%  Acquisition
%    tom_tilt_series         Tool to acquire a tilt series (WDN & FB). 
%    Note: sub-functions are not listed. See Acquisition folder or documentation
%
%  Analysis
%    tom_ccc                        calculates normalized 3d cross correlation coefficient (FF)
%    tom_cm                         calculates the center of mass (1D, 2D, 3D, AL)
%    tom_corr                       computes 3d cross correlation fuction (FF)
%    tom_create_ctf                 calculates 2D or 3D CTF (pure phase contrast) (2D, 3D, FF)
%    tom_ctf                        calculates and plots CTF (pure phase contrast) (1D, FF) 
%    tom_ctf_zero                   calculates zeros of CTF (pure phase contrast, FF)
%    tom_dev                        calculates deviation etc. for images  (2D, 3D, AF, tested FF)
%    tom_hist3d                     calculates a histogram of 3D data(SN)
%    tom_moi                        calculates the moments of inertia (3D, FF)
%    tom_orcd                       performs three-dim. orientational search (FF)
%    tom_peak                       determines the maximum/minimum value (1D, 2D, 3D, AL & GS, tested FF)
%    tom_xraycorrect                corrects X-Ray defects on EM-Images (2D, AL)
%    tom_compare                    compares 2 arrays in Fourier space (2D, 3D, MR)
%  Average 2D
%    tom_av2_particlepickergui      tool for 2D averaging of single particles (WDN & FB)
%
%  Average 3D
%
%
%  Display
%    tom_ctffit                     tool to determine the defocus (2D, WDN)
%    tom_dspcub                     visualization of a 3D volume in a gallery (3D, AF)
%    tom_embrowse                   browser for EM files (2D,3D, WDN)  
%    tom_imagesc                    similar to imagesc, transposes the image (EM style)(WDN)
%    tom_intervol                   visualization tool for 3D data loaded in memory (3D, WDN)                               
%    tom_isosurface                 tool for isosurface representation (SN)
%    tom_particles                  tool for interactive 3D particles picking (3D, WDN, tested SN)
%    tom_picker                     tool for interactive 3D particles picking loaded in memory (FF)
%    tom_volxy                      visualization tool for 3D volumes (SN)
%
%  Filtering & Transformation
%    tom_bandpass                   performs bandpass filtering of image or volume (FF)
%    tom_cut                        filtering in Fourier Space (AL)
%    tom_filter                     convolutes with spherical or quadratic kernel (FF)
%    tom_fourier                    calculates the fourier transformation of an image (1D, 2D, 3D, AL, tested FF)
%    tom_ifourier                   calculates the inverse fourier transformation of an image (1D, 2D, 3D, AL, tested FF)
%    tom_laplace                    performs Laplace filtering of image or volume (FF)
%    tom_oscar                      performs 3D cross correlation of template and volume (FF)
%    tom_ps                         calculates power spectrum of an image (2D, 3D, FF, tested FF)
%    tom_shift                      shifts image by a vector (even subpixels) (1D, 2D, 3D, FF)
%    tom_smooth                     performs smoothing of borders (2D, 3D, FF, tested FF)
%    tom_smoothz                    performs smoothing of borders along z only (2D, 3D, FF, tested FF)
%    tom_wedge                      calculates wedge - to be used as a filter for tomo (3D, FF)
%    tom_weight3d                   weighting function for tomographic reconstruction (FF & FB)
%
%  Geometrical shapes
%    tom_circle                     generates a circle (AL)
%    tom_cylinder                   creates an 3D array with a cylinder inside (SN)
%    tom_cylindermask               masks volume with (fuzzy) cylinder (3D, FF)
%    tom_ellipsemask                masks volume with (fuzzy) ellipse (3D, FF)
%    tom_sphere                     creates volume with a sphere (FF)
%    tom_spheremask                 masks volume with (fuzzy) sphere  (2D, 3D, FF)
%
%  Input/Output functions
%    tom_convert                    converts images to a different file format (AL & GS)
%    tom_emheader                   Adds a EM-header structure to a matrix (SN)
%    tom_emread                     reads data in EM-file format (SN)
%    tom_emreadc                    reads data in EM-file format. Developer version !(SN)
%    tom_emwrite                    writes data in EM-file format (SN)
%    tom_emwritec                   writes data in EM-file format. Developer version !(SN)
%    tom_isemfile                   checks if data is in EM-file format (SN)
%    tom_ismrcfile                  checks if data is in MRC-file format (SN)
%    tom_mrcreadclassic             reads data in MRC-file format (SN)
%    tom_mrcstack2emseries          converts a MRC-file format stack in a series of single EM-files (SN)
%    tom_mrcstack2emstack           converts a MRC-file format stack in a EM-format stack (SN)
%    tom_mrcwriteclassic            writes data in MRC-file format (SN)
%    tom_pdb2em                     Converts PDB-files to EM files (SN)
%    tom_pdbread                    reads a Protein Data Bank (pdb) file into a MATLAB structure (SN)
%    tom_rawread                    user interface for reading raw data (SN)
%    tom_reademheader               reads only the header of an EM-file (SN)
%    tom_readmrcheader              reads only the header of an MRC-file (SN)
%
%  Miscellaneous
%    tom_error                      generates Gaussian or Poisson noise (AL)
%    tom_makemovie                  generates an avi-file format movie file from a 3D stack (SN)
%
%  Reconstruction
%    tom_backproj3d                 performs 3D backprojection (SN)
%    tom_rec3d                      a tool for tomographic 3D reconstructions (WDN)
%    tom_reconstruction3d           performs an aligned, weighted backprojection (SN & FB)
%    tom_recparticles               to reconstruct particles with high resolution, from tom_particles (3D, WD, tested SN)
%    tom_setmark                    tool for the alignment of projections (2D, WDN, tested FF)
%
%  Spatial Transformation
%    tom_bin                        performs binning of 1D, 2D or 3D images (1D,2D,3D , FF, tested SN)
%    tom_cart2cyl                   transforms 3D-volumes from cartesian to cylinder coordinates (FF) 
%    tom_cart2polar                 transforms 2D-images from cartesian to polar coordinates (2D, FF)
%    tom_cart2sph                   transforms 3D-volumes from cartesian to spherical coordinates (3D, FF)
%    tom_cyl2cart                   transforms 3D volume from polar to cartesian coordinates(FF)
%    tom_limit                      sets limits to the values of an image (2D, 3D, AL)
%    tom_mirror                     creates the mirror of an image (2D, 3D, AL)
%    tom_move                       moves the image in x,y,z direction (2D, 3D, AL, tested FF)
%    tom_paste                      pastes array a into array b (AL)
%    tom_pointrotate                rotates a 3D vector point (FF)
%    tom_polar2cart                 transforms 2D-images from polar to cartesian coordinates (2D, FF)
%    tom_red                        extracts a subimage/subvolume from the input image/volume (AL)
%    tom_rotate                     rotates a 2D or 3D data  (3D, SN & FB)
%    tom_sph2cart                   transforms 3D-volumes from polar to cartesian coordinates (3D, FF)
%    tom_symref                     does n-fold symmetrization of a 3D reference (3D, FF)
%
%  Utilities
%    tom_chooser                    subjective classification of particles based on a motive list (3D, FF, tested FF)
%    tom_sortseries                 copy, correct x-rays and sort a tiltseries by angles (SN)
%    tom_tilt_lines                 plots marker points (FF) 
%
%  To Do:
%    tom_proj3d2d  
%    tom_phantom3d
%    tom_proj3d2d
%    tom_mrc2em
%    tom_mrcread                    reads data in MRC-file format, FEI style (SN)
%    tom_spiderread                 reads data in SPIDER-file format (SN)
%    tom_infoseries
%    tom_imadj
%    tom_symref2d                   does n-fold symmetrization of a 2D reference (2D, FF)
%    tom_tiltstrech
%    tom_unbend
%    tom_crowther                   tool for estimating the resolution with the s.c. Crowther criteria (SN)
%    tom_average3d
%    
%
%  Installation
%    makefun                        compiles Functions (please extend after creation of Functions ...)	
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
