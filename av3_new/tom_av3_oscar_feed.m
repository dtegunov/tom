function ret = tom_av3_oscar_feed(volumefile, templatefile, Outputfile, psffile, maskfile, logfile, fftsize, nr_procs, angle_start, angular_incr, angle_end, scriptflag)
%TOM_AV3_OSCAR_FEED feed oscar with all needed variables
%
%   ret = tom_av3_oscar_feed(volumefile, templatefile, Outputfile, psffile, maskfile, logfile, fftsize, nr_procs, angle_start, angular_incr, angle_end, scriptflag)
%
%PARAMETERS
%
%  INPUT
%   volumefile          full path to volume file
%   templatefile        full path to template file
%   Outputfile          full path to output files, without extension, .ccf, .ccf.norm, .ang will be added automatically
%   psffile             full path to point spread function file
%   maskfile            full path to mask file
%   logfile             full path to log file
%   fftsize             size of fft that will be processed on each cpu, must be power of two and smaller than volume dimensions
%   nr_procs            number of cpus
%   angle_start         [phi psi theta]
%   angular_incr        [phi psi theta]
%   angle_end           [phi psi theta]
%   scriptflag          ...
%  
%  OUTPUT
%   ret                 ...
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCAR, TOM_AV3_PASTE4OSCARGUI, TOM_AV3_OSCARGUI
%
%   created by AK 07/10/05
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


if nargin < 12
    scriptflag = 0;
end

ret = 1;

[a b] = unix('hostname');
if ~isempty(findstr(b,'berlin'))
   molmatch_path='/fs/pool/pool-bmsan-apps/oscar/oscar_berlin/bin/oscar';
   disp('use 64-bit version for Itanium 2');
   %!lamnodes
elseif ~isempty(findstr(b,'cluster'))
   if (~isempty(findstr(b,'09')) | ~isempty(findstr(b,'10')) | ~isempty(findstr(b,'11')) | ~isempty(findstr(b,'12')))
        molmatch_path='/raid5/apps/titan/oscar_titan/bin/oscar';
   else
        molmatch_path='/fs/pool/pool-bmsan-apps/oscar/oscar_cluster/bin/oscar';
   end;
   disp('use 64-bit version for Opteron');
   %!lamnodes
elseif ~isempty(findstr(b,'atlas')) || ~isempty(findstr(b,'prometheus')) || ~isempty(findstr(b,'calypso'))
   molmatch_path='/fs/pool/pool-bmsan-apps/oscar/oscar_titan/bin/oscar';
   disp('oscar Verson for Titan CP');
else
    %molmatch_path='/fs/pool/pool-bmsan-apps/oscar/oscar_cluster/bin/oscar';
%    molmatch_path='/fs/b_baumei/nickell/molmatch_v1.0/bin/molmatch.exe';
    molmatch_path='/fs/pool/pool-bmsan-apps/oscar/oscar_bmcluster/bin/oscar';
end;
disp(['Cpus used: ' num2str(nr_procs)]);
%molmatch_path='/fs/pool/pool-bmsan-apps/oscar/oscar_berlin/bin/oscar';

header = tom_reademheader(templatefile);
tsize = header.Header.Size;

%generate mask
if exist(maskfile,'file') == 0
    mask=ones(128,128,128);
    %mask=tom_spheremask(ones(tsize(1),tsize(2),tsize(3)),tsize(1)./2-5,5,[tsize(1)./2+1,tsize(2)./2+1,tsize(3)./2+1]);
    tom_emwrite('mask.em',mask);
    mask = 'mask.em';
    disp('No mask file found, generating default');
else
    header = tom_reademheader(maskfile);
    msize = header.Header.Size;
    mask = maskfile;
    if tsize ~= msize
        disp('Error: Template and mask dimensions must agree!');
        return;
    end;
end;

%generate psf
if exist(psffile,'file') == 0
    yyy = zeros(tsize(1),tsize(2),tsize(3));
    wedge=tom_wedge(yyy,1);
    yyy(1,1,1) =1;
    psf = real(tom_ifourier(ifftshift(fftshift(tom_fourier(yyy)).*wedge)));
    
    tom_emwrite('psf.em', psf);
    psf = 'psf.em';
    disp('No psf file found, generating default');
else
    header = tom_reademheader(psffile);
    psize = header.Header.Size;
    psf = psffile;
    if tsize ~= psize
        disp('Error: Template and psf dimensions must agree!');
        return;
    end;
end;

time_start = clock;

% run the parallel code
angles = [num2str(angle_start(1)) ' ' num2str(angle_end(1)) ' ' num2str(angular_incr(1)) ' ' num2str(angle_start(2)) ' ' num2str(angle_end(2)) ' ' num2str(angular_incr(2)) ' ' num2str(angle_start(3)) ' ' num2str(angle_end(3)) ' ' num2str(angular_incr(3))];
disp(['start oscar . . .']);time_molmatch=clock;
if scriptflag == 1
    [filename, pathname] = uiputfile('*.sh', 'Save Script as');
    if ischar(filename)
        call = ['mpirun -np ' num2str(nr_procs) ' ' molmatch_path ' ' volumefile ' ' templatefile ' ' Outputfile  ' ' angles ' ' psffile ' ' maskfile ' ' num2str(fftsize) ' | tee /dev/tty ' logfile ''];
        fid = fopen([pathname filename],'wt');
        fprintf(fid,'%s\n',call);
        fclose(fid);
        unix(['chmod +x ' pathname filename]);
    end;
else
    %[status result]=unix(['ssh berlin "cd /fs/sally02/lv02/pool/pool-nickell/thermosome/data_gutsche/new_average/auto_refinements/110407 && mpirun -ssi rpi tcp -np ' num2str(nr_procs) ' ' molmatch_path ' ' volumefile ' ' templatefile ' ' Outputfile  ' ' angles ' ' psffile ' ' maskfile ' ' num2str(fftsize) ' | tee ' logfile '"'],'-echo');
    [status result]=unix(['mpirun -ssi rpi tcp -np ' num2str(nr_procs) ' ' molmatch_path ' ' volumefile ' ' templatefile ' ' Outputfile  ' ' angles ' ' psffile ' ' maskfile ' ' num2str(fftsize) ' | tee ' logfile],'-echo');
    if status==0
        disp(['completed.']);
        try
            disp(['last oscar cycle needed: ' num2str(etime(clock,time_molmatch)./60) ' minutes' ]);
        end
    else
        disp(['problems with oscar !!!']); 
        ret = 0;
        return;
    end;
end