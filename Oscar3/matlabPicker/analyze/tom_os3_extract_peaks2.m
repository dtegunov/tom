function [coord,align2d]=tom_os3_extract_peaks2(peaks,num_of_peaks,templateSize,filename,scaling_factor,zero_rad)
%  TOM_OS3_EXTRACT_PEAKS2 extract peaks from given cc-function 
%  
%     [coord,align2d]=tom_os3_extract_peaks2(peaks,num_of_peaks,templateSize,filename,scaling_factor)
%  
%  PARAMETERS
%  
%    INPUT
%     peaks               ccf function
%     num_of_peaks        number of peaks 2 be extracted
%     templateSize        size of template must correspond 2 ccf
%     filename            filname of the org image for align2d struct
%     scaling_factor      used for upscaling the coord       
%     zero_rad            radius in percent of template/2 
%                         100 means radius=template/2  
%
%
%    OUTPUT
%     coord            extracted n coordinates (scaled)
%     align2d          align2d struct containing scaled coords
%
%  EXAMPLE
%
%  
% %200 peaks extracted from cc and scaled by a fact of 2^2
% [coord,align2d]=tom_os3_extract_peaks2(cc,300,[64 64],'/fs/pool/pool-nickell/26S/em/data/bohn/2d/090902_p47f11/high/p47f11_3680.em',2)
%
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

if (nargin < 6)
    zero_rad=100;
end;

rad=round(templateSize(1)./2);

peaks = tom_os3_eraseBorders(peaks,zeros(round(templateSize)),0);


for i=1:num_of_peaks
    [coord(i,:) value(i) peaks] = tom_peak(peaks,round(rad.*(zero_rad./100)));
end;

for pointnumber=1:num_of_peaks
    align2d(1,pointnumber).dataset = '';
    align2d(1,pointnumber).filename = filename;
    align2d(1,pointnumber).position.x = coord(pointnumber,1).*(scaling_factor);
    align2d(1,pointnumber).position.y = coord(pointnumber,2).*(scaling_factor);
    align2d(1,pointnumber).class = 'default';
    align2d(1,pointnumber).radius = round((templateSize(1).*scaling_factor)./2);
    align2d(1,pointnumber).color = [0 1 0];
    align2d(1,pointnumber).shift.x = 0;
    align2d(1,pointnumber).shift.y = 0;
    align2d(1,pointnumber).angle = 0;
    align2d(1,pointnumber).isaligned = 0;
    align2d(1,pointnumber).ccc = value(pointnumber);
    align2d(1,pointnumber).quality = 0;
    align2d(1,pointnumber).normed = 'none';
    align2d(1,pointnumber).ref_class = 0;
    coord(pointnumber,:)=coord(pointnumber,:).*scaling_factor;
end