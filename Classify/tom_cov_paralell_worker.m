function e=tom_cov_paralell_worker(infile,outfile,package,pause_time)
%TOM_COV_PARALLEL_WORKER creates ...
%
%   e=tom_cov_paralell_worker(infile,outfile,package,pause_time)
%
%PARAMETERS
%
%  INPUT
%   infile              ...
%   outfile             ...
%   package             ...
%   pause_time          ...
%  
%  OUTPUT
%   e           		...
%
%EXAMPLE
%   ... = tom_cov_paralell_worker(...);
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

error(nargchk(0, 7, nargin, 'struct'))

h=tom_reademheader(infile);
sz_m=h.Header.Size;

if nargin < 3
    package=[1 sz_m(2)];
end;

e=1;

start=package(1);
stop=package(2);

m=tom_emreadc(infile);
m=m.Value;
zz=1;

for i=package(1):package(2)
    for ii=1:(size(m,2))
        out(zz,ii)=sum( (m(:,i)-mean(m(:,i))) .*(m(:,ii)-mean(m(:,ii))) ) ./ (size(m(:,1),1)-1);

    end;
    zz=zz+1;
end;

pause(pause_time);
tom_emwritec(outfile,single(out),'subregion',[(package(1)) 1 1],[size(out,1) size(out,2) 1]);


