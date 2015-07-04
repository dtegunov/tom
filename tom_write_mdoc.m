function [ ] = tom_write_mdoc( mdoc, filename )

fid = fopen(filename, 'W');

fprintf(fid, 'PixelSpacing = %f\n', mdoc.pixelsize);
fprintf(fid, 'ImageFile = %s\n', mdoc.imagefile);
fprintf(fid, 'ImageSize = %d %d\n', mdoc.size(1), mdoc.size(2));
fprintf(fid, 'DataMode = 1\n');
fprintf(fid, '\n');
fprintf(fid, '[T = SerialEM]\n');
fprintf(fid, '\n');
fprintf(fid, '[T = Tilt axis angle = %f, binning = %d]\n', mdoc.axisangle, mdoc.binning);
fprintf(fid, '\n');

for i = 1:size(mdoc.zvalue, 2)
    fprintf(fid, '[ZValue = %d]\n', mdoc.zvalue(i));
    if isfield(mdoc, 'theta')
        fprintf(fid, 'TiltAngle = %f\n', mdoc.theta(i));
    end;
    if isfield(mdoc, 'defocus')
        fprintf(fid, 'Defocus = %f\n', mdoc.defocus(i));
    end;
    if isfield(mdoc, 'rot')
        fprintf(fid, 'RotationAngle = %f\n', mdoc.rot(i));
    end;
    if isfield(mdoc, 'shift')
        fprintf(fid, 'ImageShift = %f %f\n', mdoc.shift(1,i), mdoc.shift(2,i));
    end;
    if isfield(mdoc, 'exptime')
        fprintf(fid, 'ExposureTime = %f\n', mdoc.exptime(i));
    end;
    if isfield(mdoc, 'expdose')
        fprintf(fid, 'ExposureDose = %f\n', mdoc.expdose(i));
    end;
    if isfield(mdoc, 'cummexpdose')
        fprintf(fid, 'CummulativeExposureDose = %f\n', mdoc.cummexpdose(i));
    end;
    if isfield(mdoc, 'numframes')    
        fprintf(fid, 'NumSubFrames = %d\n', mdoc.numframes(i));
    end;
    if isfield(mdoc, 'timestamp')    
        fprintf(fid, 'DateTime = %s\n', datestr(mdoc.timestamp(i), 'dd-mmm-yy HH:MM:SS'));
    end;
    fprintf(fid, '\n');
end;

fclose(fid);

end

