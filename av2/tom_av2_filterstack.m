function tom_av2_filterstack(instackfile,inalignfile,field,operator,value,outstackfile,outalignfile)
%TOM_AV2_FILTERSTACK filters a particle stack
%
%   tom_av2_filterstack(instackfile,inalignfile,field,operator,value,outstackfile,outalignfile)
%
%PARAMETERS
%
%  INPUT
%   instackfile         full path to input stack file (.em)
%   inalignfile         full path to input alignment file (.mat)
%   field               fieldname in align2d structure to filter
%   operator            relational operator (eg. ==)
%   value               value to compare the field with
%   outstackfile        full path to output stack file (.em)
%   outalignfile        full path to output alignment file (.mat)
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_filterstack(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER
%
%   created by AK 01/13/06
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


error(nargchk(7,7,nargin));

if ~isstruct(inalignfile)
    try
        s = load(inalignfile);
    catch
        error('Could not load alignment file.');
    end
else
    s.align2d = inalignfile;
end

if tom_isemfile(instackfile) ~= 1
    error('Could not load stack file.');
end

h = tom_reademheader(instackfile);
newstack = [];
c = 1;

stack_sz=0;
%check out size of new stack
for i=1:size(s.align2d,2)
    t = eval(['s.align2d(i).(field) ' operator ' value;']);
    if t == 1
        stack_sz=stack_sz+1;
    end;
end;

if stack_sz > 0
    %allocate some memory
    newstack=zeros(h.Header.Size(1),h.Header.Size(2),stack_sz);

    zz=1;
    for i=1:size(s.align2d,2)
        t = eval(['s.align2d(i).(field) ' operator ' value;']);
        if t == 1
            align2d(c) = s.align2d(i);
            c = c+1;
            im = tom_emreadc(instackfile,'subregion',[1 1 i],[h.Header.Size(1)-1 h.Header.Size(2)-1 0]);
            %newstack = cat(3,newstack,im.Value);
            newstack(:,:,zz)=im.Value;
            zz=zz+1;
        end
    end
    try;align2d = rmfield(align2d,'selectedtowrite');end;
    save(outalignfile,'align2d');
    tom_emwrite(outstackfile,newstack);
else
    error('New stack would contain no particles, please check your input.');
end