% Used to convert an bits and byts from the em-header to an em-structure 
% (with the fields 'Value' and 'Header').
function st = tom_emheader_data_2_struct(data, em_size, magic, comment, emdata, userdata, filename)


Header = struct();
Header.Size = size(data);
Header.Voltage = double(emdata(1));
Header.Cs = double(emdata(2)) / 1000;
Header.Aperture = double(emdata(3));
Header.Magnification = double(emdata(4));
Header.Postmagnification = double(emdata(5)) / 1000;
Header.Exposuretime = double(emdata(6)) / 1000;
Header.Objectpixelsize = double(emdata(7)) / 1000;
Header.Microscope = tom_em_microscope_name(emdata(8));
Header.Pixelsize = double(emdata(9)) / 1000;
Header.CCDArea = double(emdata(10)) / 1000;
Header.Defocus = double(emdata(11));
Header.Astigmatism = double(emdata(12));
Header.AstigmatismAngle = double(emdata(13)) / 1000;
Header.FocusIncrement = double(emdata(14));
Header.CountsPerElectron = double(emdata(15)) / 1000;
Header.Intensity = double(emdata(16)) / 1000;
Header.EnergySlitwidth = double(emdata(17));
Header.EnergyOffset = double(emdata(18));
Header.Tiltangle = double(emdata(19)) / 1000;
Header.Tiltaxis = double(emdata(20)) / 1000;
Header.Marker_X = double(emdata(23));
Header.Marker_Y = double(emdata(24));
Header.Username = char(userdata(1:20));
Header.Date = char(userdata(21:28));
Header.Comment = comment;

warning off

if (exist('filename','var') && ~isempty(filename))
    [Header.Pathname, Header.Filename, ext] = fileparts(getAbsFilename(filename));
    Header.Filename = [Header.Filename, ext];
else 
    Header.Pathname = '';
    Header.Filename = '';
end;

warning on;

Header.EM.Size = em_size;
Header.EM.Magic = magic;
Header.EM.Comment = comment;
Header.EM.Parameter = emdata;
Header.EM.Fillup = userdata;

st = struct('Value', {data}, 'Header', {Header});




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = getAbsFilename(s)


[success, message] = fileattrib(s);
if (success && ~message.directory)
    s = message.Name;
end;




