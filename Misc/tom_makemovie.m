function tom_makemovie(in, out, fps, quality, range, average)
%TOM_MAKEMOVIE generates an avi-file format movie file from a 3D array
%
%   tom_makemovie(in, out, fps, quality, range, average)
%
%   Scroll in z-direction and write out movie from single slices or 
%   average along the z-axis. Nice for visualizing 3D volumes of
%   a tomographic volume.
%
%PARAMETERS
%
%  INPUT
%   in                  3D array
%   out                 Filename of the generated output file
%   fps                 frames per second in the movie
%   quality             compression rate (100% without loss)
%   range               scaling range of the data in a vector
%   average             average along the z-axis multiple layers
%  
%  OUTPUT
%   writes a movie directly to the harddrive
%
%EXAMPLE
%   tom_makemovie(in,'Movie.avi',10,100,[-1.5 1.5],3);
%
%   Make a movie 'Movie.avi' with 10 frames per second and 100 per cent
%   quality. Use -1.5 as a lower limit for the gray values and 1.5
%   as an upper limit and average 3 layers in z-direction
%
% REMARK
%   
%   Check the different compression schemes in line 48:  
%   'Cinepak', 'Indeo3', 'Indeo5', 'MSVC',', 'RLE', 'None'
%   maybe there are different ones installed at your Windows
%   distribution
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 01/03/03
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


fig=figure;set(gcf,'Position',[50 50 850 850]);

set(fig,'DoubleBuffer','on');
mov = avifile(out);
mov.Fps=fps;
mov.Compression='None';
% 'Cinepak', 'Indeo3', 'Indeo5', 'MSVC',', 'RLE', 'None'
mov.Quality=quality;

for lauf=120:210-average; 
    im=mean(in(:,:,lauf:lauf+average),3);
    imagesc(im',[range(1) range(2)]);
    axis image;
    axis off;
    colormap gray;
    F = getframe(gca);
    mov = addframe(mov,F);
end
mov = close(mov);
close(fig);