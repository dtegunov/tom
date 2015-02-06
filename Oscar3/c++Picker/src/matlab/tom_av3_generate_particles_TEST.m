%!rm -f /fs/home/haller/haller/DA/data/gen_part/* /fs/home/haller/haller/DA/data/outputdir/*     
outputdir = '/fs/home/haller/haller/DA/data/gen_part';
orig_template = '/fs/home/haller/haller/DA/data/Ribo_0_68_nm.em';
wedge_template = [];
wedge_particles = [30,25];
nparticles = 1;
do_rotate = true;
sigma_shift = 10;
round_shift = true;
sigma_noise = 0;
append_angles_list = 0;


ss = size(getfield(tom_emreadc3(orig_template), 'Value'));
if (1)
    if (~isempty(wedge_particles))
        wedge = tom_wedge(ones(ss), wedge_particles(1));
        if (wedge_particles(2) > 0) 
            wedge = wedge .* tom_spheremask(ones(ss), wedge_particles(2));
        end;
        wedge_particles = wedge;
    end;
end;

if (1)
    if (~isempty(wedge_template))
        wedge = tom_wedge(ones(ss), wedge_template(1));
        if (wedge_template(2) > 0) 
            wedge = wedge .* tom_spheremask(ones(ss), wedge_template(2));
        end;
        wedge_template = wedge;
    end;
end;

 
rand('state', 2);

tom_av3_generate_particles(outputdir, orig_template, wedge_template, wedge_particles, nparticles, do_rotate, sigma_shift, round_shift, sigma_noise, append_angles_list);
!ssh calypso make -C ~/haller/DA/tomc x_gen



[a, binning] = system('sed -n ''s/^binning=\(.*\)$/\1/p'' /fs/home/haller/haller/DA/tomc/data/corr/config.txt'); 
if (a)
    error('system');
end;
clear a;
binning = str2double(binning);
if (round(binning)~=binning || ~isfinite(binning) || binning<1)
    binning = 1;
end;


template_original = tom_emreadc3(fullfile(outputdir, 'template_original.em'), [], [], binning*[1 1 1]); template_original = template_original.Value;
ss = size(template_original);

particle_wedge = tom_emreadc3(fullfile(outputdir, 'particle_wedge.em'), [], [], binning*[1 1 1]); particle_wedge = particle_wedge.Value;
template_wedge = tom_emreadc3(fullfile(outputdir, 'template_wedge.em'), [], [], binning*[1 1 1]); template_wedge = template_wedge.Value;

particle = tom_emreadc3(fullfile(outputdir, 'particle_001.em'), [], [], binning*[1 1 1]); particle = particle.Value;
template = tom_emreadc3(fullfile(outputdir, 'template.em'), [], [], binning*[1 1 1]); template = template.Value;

angles = tom_emreadc3(fullfile(outputdir, 'angles.em')); angles = angles.Value * 180 / pi;
shifts = tom_emreadc3(fullfile(outputdir, 'shifts.em')); shifts = shifts.Value;


if (~exist('ccv','var') || ndims(ccv)~=ndims(ccv_old) || any(ss~=size(ccv_old)))
    ccv = nan(ss);
end;
ccv_old = ccv;


outputdir2 = '/fs/sandy01/lv03/pool/bmsan/haller/DA/data/outputdir';
ccv = tom_emreadc3(fullfile(outputdir2, 'ccv_____0_____0.em')); ccv = fftshift(ccv.Value);


if (~exist('ccvval','var'))
    ccvval = nan;
end;
ccvval_old = ccvval;


[ccvpeak, ccvval] = tom_peak(ccv); ccvpeak = ccvpeak - (floor(size(ccv,1)/2)+1);

disp(['peak: [' num2str(ccvpeak) ']: ' num2str(shifts) '      ' num2str(ccvval)]);
disp(num2str([ccvval, ccvval_old, ccvval - ccvval_old]));
disp(num2str(angles'));
tom_dev(abs(ccv - ccv_old));

if (0)
a1 = tom_emreadc3('/fs/home/haller/haller/DA/data/outputdir/WWwedge.em');  a1 = a1.Value;
a2 = tom_emreadc3('/fs/home/haller/haller/DA/data/outputdir/WWwedge2.em'); a2 = a2.Value;
a3 = tom_emreadc3('/fs/home/haller/haller/DA/data/outputdir/WWwedge3.em'); a3 = a3.Value;
end


!cat /fs/home/haller/haller/DA/data//outputdir/peakfile.txt
disp(tom_peak(ifftshift(ccv))-1);


