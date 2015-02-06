function Align = tom_av3_extract_oscar_anglesshifts(oscarvol, radius, number_of_peaks, waitbarflag)
%TOM_AV3_EXTRACT_OSCAR_ANGLESSHIFTS extracts the angles and shifts from an oscar output volume
%
%   Align = tom_av3_extract_oscar_anglesshifts(oscarvol, radius, number_of_peaks, waitbarflag)
%
%PARAMETERS
%
%  INPUT
%   oscarvol            full path to oscar volume without extension
%   radius              radius for peak detection, typically particle size / 2
%   number_of_peaks     number of extracted peaks
%   waitbarflag         1 = GUI waitbar
%                       0 = text messages
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   Align = tom_av3_extract_oscar_anglesshifts('./cyl_search_ccf', 16, 500, 1)
%    extracts 500 peaks from the oscar CCF volume ./cyl_search_ccf.ccf.norm,
%    using a 16 pixel radius for peak detection and shows the waitbar
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCARGUI, TOM_AV3_OSCARGUI, TOM_AV_EXTRACT_ANGLESSHIFTS,
%   TOM_AV_EXTRACTANGLESSHIFTSGUI, TOM_AV3_AVERAGE
%
%   created by SN 09/20/06
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


 if nargin < 4
     waitbarflag = 0;
 end
% 
% try
%     tsize = Align(1,1).Tomogram.Header.Size';
% catch
%     prompt = {'Enter templatesize x:','Enter templatesize y:','Enter templatesize z:'};
%     dlg_title = 'Missing template size';
%     num_lines = 1;
%     def = {'','',''};
%     answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
%     tsize = [str2num(answer{1}), str2num(answer{2}), str2num(answer{3})];
% end
% 
 if waitbarflag == 1
     h = waitbar(0,'Extract peak list');
 end


header = tom_reademheader(strcat(oscarvol,'.ccf.norm'));
vsize = header.Header.Size;

%Extract angles from Comment
angs = sscanf(char(header.Header.Comment'),'Angular range: %f %f %f %f %f %f %f %f %f');
angle_start = [angs(1) angs(4) angs(7)];
angle_end = [angs(2) angs(5) angs(8)];
angular_incr = [angs(3) angs(6) angs(9)];

meanccc = [];
ccf = tom_emreadc(strcat(oscarvol,'.ccf.norm'));

for i = 1:number_of_peaks

    [peakposition peakval ccf.Value] = tom_peakc(ccf.Value,radius);
    
    ang = tom_emreadc(strcat(oscarvol,'.ang'), 'subregion',[peakposition(1) peakposition(2) peakposition(3)],[0 0 0]);

    Align(1,i).Oscar.CCF_volume = strcat(oscarvol,'.ccf.norm'); 
    Align(1,i).Tomogram = strcat(oscarvol,'.ccf.norm'); 
    Align(1,i).Oscar.Angle_volume = strcat(oscarvol,'.ang'); 
    Align(1,i).Oscar.Angular_range_start = [angs(1) angs(4) angs(7)]; 
    Align(1,i).Oscar.Angular_range_stop = [angs(2) angs(5) angs(8)]; 
    Align(1,i).Oscar.Angular_range_step = [angs(3) angs(6) angs(9)]; 
    Align(1,i).Oscar.Position.X = (peakposition(1)); 
    Align(1,i).Oscar.Position.Y = (peakposition(2));
    Align(1,i).Oscar.Position.Z = (peakposition(3));
    Align(1,i).Oscar.CCC = peakval;    
    Align(1,i).CCC = peakval;    
    meanccc = [meanccc peakval];
    
    angval = ang.Value;
    [angle_out] = av3_index2angle_new(angval, angle_start, angular_incr, angle_end);
    Align(1,i).Oscar.Angle.Phi = -angle_out(2);
    Align(1,i).Oscar.Angle.Psi = -angle_out(1);
    Align(1,i).Oscar.Angle.Theta = -angle_out(3);
    Align(1,i).Oscar.Angle.Rotmatrix = tom_angles2rotmatrix([-angle_out(2) -angle_out(1) -angle_out(3)]);
    Align(1,i).Angle.Phi = -angle_out(2);
    Align(1,i).Angle.Psi = -angle_out(1);
    Align(1,i).Angle.Theta = -angle_out(3);
    Align(1,i).Angle.Rotmatrix = tom_angles2rotmatrix([-angle_out(2) -angle_out(1) -angle_out(3)]);
    Align(1,i).Shift.X=0;
    Align(1,i).Shift.Y=0;
    Align(1,i).Shift.Z=0;
    Align(1,i).Class=0;
    Align(1,i).NormFlag=0;
    Align(1).Filter=[0 0];
    Align(1).ProjectionClass=0;
    
    if waitbarflag == 1
        waitbar(i./number_of_peaks,h,[num2str(i), ' of ', num2str(number_of_peaks), ' particles done.']);
    elseif waitbarflag == 0
        disp([num2str(i), ' of ', num2str(number_of_peaks), ' particles done.']);
    end
end



if waitbarflag == 1
    close(h);
end

%do some sanity checking

if mean(meanccc,2) < 0.15
   if waitbarflag == 1
        warndlg('WARNING: mean CCC of this run is below 0.15, the quality of this oscar run may be very low!');
   else
       disp('WARNING: mean CCC of this run is below 0.15, the quality of this oscar run may be very low!');
   end
end