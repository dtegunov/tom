function [ mdoc ] = tom_read_mdoc( filename )

fid = fopen(filename);

header = textscan(fid, '%s', 8, 'Delimiter', '\n');
header = header{1};
mdoc.pixelsize = equalval(header{1}, '%f');
mdoc.imagefile = equalval(header{2}, '%s');
mdoc.size = equalval(header{3}, '%f %f');
mdoc.axisangle = valinline(header{8}, 'Tilt axis angle', '%f');
mdoc.binning = valinline(header{8}, 'binning', '%f');

mdoc.zvalue = [];
mdoc.numframes = [];
mdoc.rot = [];
mdoc.theta = [];
mdoc.defocus = [];
mdoc.exptime = [];
mdoc.expdose = [];
mdoc.cummexpdose = [];
mdoc.shift = [];
mdoc.timestamp = [];

while (~feof(fid))
    line = fgetl(fid);
    while (~feof(fid) && size(strfind(line, 'ZValue'), 1) == 0)
        line = fgetl(fid);
    end;
    
    if size(strfind(line, 'ZValue'), 1) == 0
        break;
    end;
    
    mdoc.zvalue = [mdoc.zvalue, sscanf(line, '[ZValue = %f]')];
    
    while (~isempty(line))
        parts = strsplit(line, ' = ');
        if strcmp(parts{1}, 'TiltAngle')
            mdoc.theta = [mdoc.theta, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'Defocus')
            mdoc.defocus = [mdoc.defocus, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'RotationAngle')
            mdoc.rot = [mdoc.rot, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'ImageShift')
            mdoc.shift = [mdoc.shift, sscanf(parts{2}, '%f %f')];
        elseif strcmp(parts{1}, 'ExposureTime')
            mdoc.exptime = [mdoc.exptime, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'ExposureDose')
            mdoc.expdose = [mdoc.expdose, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'CummulativeExposureDose')
            mdoc.cummexpdose = [mdoc.cummexpdose, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'NumSubFrames')
            mdoc.numframes = [mdoc.numframes, sscanf(parts{2}, '%f')];
        elseif strcmp(parts{1}, 'DateTime')
            datestring = sscanf(parts{2}, '%s');
            mdoc.timestamp = [mdoc.timestamp, datenum(datestring, 'dd-mmm-yyHH:MM:SS')];
        end;
        
        if feof(fid)
            break;
        end;
        line = fgetl(fid);
    end;
end;
fclose(fid);

names = fieldnames(mdoc);
for n=1:numel(names)
    if isempty(mdoc.(names{n}))
        mdoc = rmfield(mdoc, names{n});
    end;
end;

end

function v = equalval ( text, valformat )
    parts = strsplit(text, ' = ');
    v = parts{2};
    v = sscanf(v, valformat);
end

function v = valinline ( line, valname, valformat )
    pos = strfind(line, valname);
    line = line(pos:end);
    v = equalval(line, valformat);
end