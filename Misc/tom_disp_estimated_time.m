function [store]=tom_disp_estimated_time(store,flag,name)
%TOM_DISP_ESTIMATED_TIME creates ...
%
%   [store]=tom_disp_estimated_time(store,flag,name)
%
%PARAMETERS
%
%  INPUT
%   store               ...
%   flag                ...
%   name                ...
%  
%  OUTPUT
%   store               ...
%
%EXAMPLE
%   ... = tom_disp_estimated_time(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

% funny pseudo constructor
if (strcmp(flag,'start')) 
    fprintf('\n%s:\n  ',name);
    store.t0 = clock;
    store.t_run = clock;
    store.i=1;
    store.i_alt=1;
    store.nn=1;
end;


if (strcmp(flag,'estimate_time'))
    if (store.i>=store.num_of_mesure & store.mesured==0)
        run1_time=round(etime(clock,store.t_run));
        restrun_time=datenum(((store.end-store.i)).*(run1_time./store.i));
        fprintf('%s\n',['Expected end time: ' datestr(datevec(now)+datevec(restrun_time.*0.00001.*1./.864))]);
        drawnow;
        store.t_run = clock;
        store.mesured=1;
        return;
    end;
end;

if (strcmp(flag,'progress')) 
    if (store.i > store.i_alt)
        fprintf('%s','.')

        if ((store.i./(100.*store.nn)) >=1)
            fprintf('%d \n',store.i);
            store.nn=store.nn+1;
        end;

    end;
    
    if (store.i==store.end)
        fprintf('%s \n','end');
    end;
    store.i_alt=store.i;
end;
