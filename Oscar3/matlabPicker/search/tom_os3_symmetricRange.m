%%
function res = tom_os3_symmetricRange(img,flags)
%determine the symmetry of template and return angle range for scan

if(isstruct(img))
    img = img.Value;
end;

dim = tom_os3_fftSelect(img);
res = flags;

if(dim == 2)
    
    res.angle.start = 0;
    res.angle.end = 179;
    
    if(strcmp(flags.mode.type,'demo'))
        res.angle.end = 360;
    end;
elseif(dim == 3)
    
%     res.angle.start = 0;
%     res.angle.end = 179;    
    
    
    if(strcmp(flags.mode.type,'demo'))
        res.angle.end = 360;
    end;
    
end;