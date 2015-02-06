function tom_IMG_series2avi(basepath,ext,avi_path)
% tom_IMG_series2avi creates a movie from a image series
%  
%     tom_IMG_series2avi(basepath,avi_path,ext)
%  
%  PARAMETERS
%  
%    INPUT
%     basepath            basepath of the images
%     ext                 extension 
%     im                  filename of the movie
%
%  
%  EXAMPLE
%
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by SN/FB 01/24/06
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

dd=dir([basepath '*' ext]);

aviobj = avifile(avi_path);
fig=figure;

aviobj=set(aviobj,'Fps',7);

for ii=1:1
for i=1:length(dd)
    im=imread([basepath num2str(i) ext]);
    imagesc(im);
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
end;
end;
close(fig)
aviobj = close(aviobj);

