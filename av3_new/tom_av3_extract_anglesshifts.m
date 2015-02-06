function Align = tom_av3_extract_anglesshifts(oscarvol, Align, mask, waitbarflag)
%TOM_AV3_EXTRACT_ANGLESSHIFTS extract the angles and shifts from an oscar output volume
%
%   Align = tom_av3_extract_anglesshifts(oscarvol, Align, mask, waitbarflag)
%
%PARAMETERS
%
%  INPUT
%   oscarvol            full path to oscar volume without extension
%   Align               particle alignment structure
%   mask                mask to use for finding the peak values
%   waitbarflag         1 = GUI waitbar
%                       0 = text messages
%  
%  OUTPUT
%   Align               ...
%
%EXAMPLE
%   ... = tom_av3_extract_anglesshifts(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCARGUI, TOM_AV3_OSCARGUI, TOM_AV_EXTRACT_ANGLESSHIFTS,
%   TOM_AV_EXTRACTANGLESSHIFTSGUI, TOM_AV3_AVERAGE
%
%   created by AK 11/10/05
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

try
    tsize = Align(1,1).Tomogram.Header.Size';
catch
    prompt = {'Enter templatesize x:','Enter templatesize y:','Enter templatesize z:'};
    dlg_title = 'Missing template size';
    num_lines = 1;
    def = {'','',''};
    answer = inputdlg(prompt,dlg_title,num_lines,def,'on');
    tsize = [str2num(answer{1}), str2num(answer{2}), str2num(answer{3})];
end

if waitbarflag == 1
    h = waitbar(0,'Updating alignment list');
end


header = tom_reademheader(strcat(oscarvol,'.ccf.norm'));
vsize = header.Header.Size;

%Extract angles from Comment
angs = sscanf(char(header.Header.Comment'),'Angular range: %f %f %f %f %f %f %f %f %f');
angle_start = [angs(1) angs(4) angs(7)];
angle_end = [angs(2) angs(5) angs(8)];
angular_incr = [angs(3) angs(6) angs(9)];

%Calculate template dimensions
xbuffer = floor((vsize(1) ./ tsize(1)) - 1);
ybuffer = floor((vsize(2) ./ tsize(2)) - 1);
zbuffer = floor((vsize(3) ./ tsize(3)) - 1);

borderx = floor((vsize(1)-xbuffer*tsize(1))./2);
bordery = floor((vsize(2)-xbuffer*tsize(2))./2);
borderz = floor((vsize(3)-zbuffer*tsize(3))./2);

x = borderx; y = bordery; z = borderz;

xcounter = 1;
ycounter = 1;

numimages = xbuffer .* ybuffer .* zbuffer;

run = size(Align,1);
meanccc = [];

for i = 1:size(Align,2)

    ccf = tom_emreadc(strcat(oscarvol,'.ccf.norm'), 'subregion',[x+1 y+1 z+1],[tsize(1)-1, tsize(2)-1, tsize(3)-1]);
    try
        ccf.Value = ccf.Value .* mask;
    catch
        if waitbarflag == 1
        errordlg('Dimensions of template and mask file differ!','Abort');
        return;
        else
            error('Dimensions of template and mask file differ!');
        end
    end
    ang = tom_emreadc(strcat(oscarvol,'.ang'), 'subregion',[x+1 y+1 z+1],[tsize(1), tsize(2), tsize(3)]);

    [peakposition peakval] = tom_peak(ccf.Value,'spline');
    Align(run,i).Shift.X = -(peakposition(1) - tsize(1)./2 - 1); % changed !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Align(run,i).Shift.Y = -(peakposition(2) - tsize(2)./2 - 1);
    Align(run,i).Shift.Z = -(peakposition(3) - tsize(3)./2 - 1);
    Align(run,i).CCC = peakval;
    peakposition=round(peakposition);
    meanccc = [meanccc peakval];
    
    angval = ang.Value(peakposition(1),peakposition(2),peakposition(3));
    [angle_out] = tom_av3_index2angle_new(angval, angle_start, angular_incr, angle_end);
    Align(run,i).Angle.Phi = -angle_out(2);
    Align(run,i).Angle.Psi = -angle_out(1);
    Align(run,i).Angle.Theta = -angle_out(3);
      
    if xcounter < xbuffer
        x = x + tsize(1);
        xcounter = xcounter + 1;
    else
        x = borderx;
        y = y + tsize(2);
        xcounter = 1;
        ycounter = ycounter + 1;
        if ycounter > xbuffer
            ycounter = 1;
            y = bordery;
            z = z + tsize(3);
        end
    end

    if waitbarflag == 1
        waitbar(i./numimages,h,[num2str(i), ' of ', num2str(numimages), ' particles done.']);
    elseif waitbarflag == 0
        disp([num2str(i), ' of ', num2str(numimages), ' particles done.']);
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