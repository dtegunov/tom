function tom_av3_generate_particles(outputdir, orig_template, wedge_template, wedge_particles, nparticles, do_rotate, sigma_shift, round_shift, sigma_noise, append_angles_list, mask_for_snr)


if (~exist(outputdir,'dir'))
    error(['The directory "' outputdir '" does not exist']);
end;


flog = fopen(fullfile(outputdir, 'log.txt'), 'w');
fprintf(flog, '# Generating %d particles (%s)\n', nparticles, datestr(now));


% Load the template.
if (ischar(orig_template) && exist(orig_template,'file'))
    fprintf(flog, '# original template from file "%s"\n', orig_template);
    orig_template = tom_emreadc3(orig_template);
    orig_template = double(orig_template.Value);
elseif (isstruct(orig_template) && isfield(orig_template, 'Value'))
    orig_template = double(orig_template.Value);
    fprintf(flog, '# template as volume\n');
end;
if (~isnumeric(orig_template))
    error('Wrong template given (either em-filename, em-structure or 3d-matrix)');
end;
fprintf(flog, '# Size = [%dx%dx%d]"\n', size(orig_template,1), size(orig_template,2), size(orig_template,3));
fprintf(flog, '# Values: [%20.15f ~ %20.15f ~ %20.15f]\n# with standard-deviation %20.15f (sample standard deviation %20.15f)\n', min(orig_template(:)), mean(orig_template(:)), max(orig_template(:)), std(orig_template(:), 1), std(orig_template(:),0));






% Apply the wedge to the template and save it.
f = fopen(fullfile(outputdir, 'templates_list.txt'), 'w');
filename = fullfile(outputdir, 'template.em');
filename_2 = fullfile(outputdir, 'template_wedge.em');
fprintf(f, '%s\n', filename);
if (isempty(wedge_template))
    orig_template_wedged = orig_template;
    wedge = ones(size(orig_template));
    fprintf(f, 'nowedge\n');
    fprintf(flog, '# use no wedge\n');
elseif (isnumeric(wedge_template) && numel(wedge_template)<=2)
    wedge = tom_wedge(ones(size(orig_template)), wedge_template(1));
    if (numel(wedge_template)==2) 
        wedge = wedge .* tom_spheremask(ones(size(orig_template)), wedge_template(2));
    else
        wedge_template(2) = 0;
    end;    
    fprintf(flog, '# use wedge: simple %15.10f %15.10f\n', wedge_template(1), wedge_template(2));
    orig_template_wedged = ifftn(fftn(orig_template).*ifftshift(wedge));    
    fprintf(f, 'simple %15.10f %15.10f\n', wedge_template(1), wedge_template(2));
elseif (ndims(wedge_template) == ndims(orig_template) && all(size(wedge_template)==size(orig_template)) && (isnumeric(wedge_template) || islogical(wedge_template)))
    wedge = double(wedge_template);
    orig_template_wedged = ifftn(fftn(orig_template).*ifftshift(wedge));    
    fprintf(f, 'emfile %s\n', filename_2);
    fprintf(flog, '# use wedge from volume\n');
elseif (ischar(wedge_template) && exist(wedge_template, 'file'))
    fprintf(flog, '# use wedge from file %s\n', wedge_template);
    wedge = tom_emreadc3(wedge_template); wedge = double(wedge.Value);
    if (ndims(wedge) ~= ndims(orig_template) || all(size(wedge)~=size(orig_template)))
        error('The wedge from file has the wrong size.');
    end;
    orig_template_wedged = ifftn(fftn(orig_template).*ifftshift(wedged));    
    fprintf(f, 'emfile %s\n', filename_2);
else
    error('wrong wedge.');
end;
fclose(f);
fprintf(flog, '# save copy of original template to %s\n', fullfile(outputdir, 'template_original.em'));
tom_emwritec3(fullfile(outputdir, 'template_original.em'), orig_template);
fprintf(flog, '# save template with wedge to %s\n', filename);
tom_emwritec3(filename, orig_template_wedged);
fprintf(flog, '# save the wedge of the template to %s\n', filename_2);
tom_emwritec3(filename_2, wedge);



if (~exist('mask_for_snr', 'var'))
    mask_for_snr = [];
end;
    



% Create the wedge of the particle...
filename_2 = fullfile(outputdir, 'particle_wedge.em');
if (isempty(wedge_particles))
    wedge = ones(size(orig_template));
    wedge_s = 'nowedge';
    fprintf(flog, '# no wedge for particle\n');    
elseif (isnumeric(wedge_particles) && numel(wedge_particles)<=2)
    wedge = tom_wedge(ones(size(orig_template)), wedge_particles(1));
    if (numel(wedge_particles)==2) 
        wedge = wedge .* tom_spheremask(ones(size(orig_template)), wedge_particles(2));
    else
        wedge_particles(2) = 0;
    end;    
    wedge_s = sprintf('simple %15.10f %15.10f', wedge_particles(1), wedge_particles(2));
    fprintf(flog, '# wedge for particle: simple %15.10f %15.10f\n', wedge_particles(1), wedge_particles(2));
elseif (ndims(wedge_particles) == ndims(orig_template) && all(size(wedge_particles)==size(orig_template)) && (isnumeric(wedge_particles) || islogical(wedge_particles)))
    wedge = double(wedge_particles);
    wedge_s = sprintf('emfile %s', filename_2);
    fprintf(flog, '# wedge for particle as volume\n');    
elseif (ischar(wedge_particles) && exist(wedge_particles, 'file'))
    fprintf(flog, '# wedge for particle from emfile %s\n', wedge_particle);    
    wedge = tom_emreadc3(wedge_particles); wedge = double(wedge.Value);
    if (ndims(wedge) ~= ndims(orig_template) || all(size(wedge)~=size(orig_template)))
        error('The wedge from file has the wrong size.');
    end;
    wedge_s = sprintf('emfile %s', filename_2);
else
    error('wrong wedge.');
end;
fprintf(flog, '# save the wedge of the particle to %s\n', filename_2);
tom_emwritec3(filename_2, wedge);



% write log file.
if (~exist('do_rotate', 'var') || isempty(do_rotate) || ~do_rotate)
    fprintf(flog, '# do_rotate: no\n');
    do_rotate = false;
else
    fprintf(flog, '# do_rotate: yes\n');
    do_rotate = true;
end;
if (~exist('sigma_shift', 'var') || sigma_shift<=0)
    sigma_shift = 0;
    fprintf(flog, '# no shift\n');
else
    fprintf(flog, '# shift: sigma = %20.15f\n', sigma_shift);
end;
if (~exist('round_shift', 'var') || isempty(round_shift) || ~round_shift)
    fprintf(flog, '# round_shift: no\n');
    round_shift = false;
else
    fprintf(flog, '# round_shift: yes\n');
    round_shift = true;
end;
if (~exist('sigma_noise', 'var') || sigma_noise<=0)
    sigma_noise = 0;
    fprintf(flog, '# Add noise: no\n');
else
    fprintf(flog, '# Add noise: yes (stddev=%20.15f)\n', sigma_noise);
end;
if (isempty(mask_for_snr))
    fprintf(flog, '# don''t use mask for snr\n');
else
    ffff = fullfile(outputdir, 'mask_for_snr.em');
    fprintf(flog, ['# use mask for snr: save to "' ffff '"\n'], ffff);
    mask_for_snr = mask_for_snr ~= 0;
    tom_emwritec3(ffff, mask_for_snr);
end;
fprintf(flog, '# Shift has a radius with sigma of %20.15f\n', sigma_shift);
fprintf(flog, '# Output directory: "%s"\n\n', outputdir);





% generate the particles...
partlistname = fullfile(outputdir, 'particles_list.txt');
fp = fopen(partlistname, 'w');

angles = nan(nparticles, 3);
shifts = nan(nparticles, 3);

ss = size(orig_template);

mtemplate = mean(orig_template(:));
stemplate = std(orig_template(:), 1);
%template = (orig_template - mtemplate) / stemplate;
template = orig_template;

is_wedge_particle = any(wedge(:) ~= 1);

wedge = ifftshift(wedge);


for (i=1:nparticles)
    
    particle = template;
    
    angle = [0,0,0];
    if (do_rotate)
        angle = tom_rotmatrix2angles(rand_rot());
        particle = tom_rotatec2(particle, angle);
    end;
    angles(i, 1:3) = angle;    
    
    shift = [0,0,0];
    if (sigma_shift > 0)
        shift = rand_rot() * randn(3,1) * sigma_shift;
        if (round_shift) 
            shift = round(shift);
        end; 

        particle = tom_shift(particle, shift);
        idx = true(ss);
        idx(max(1,1+round(shift(1))):min(ss(1),ss(1)+round(shift(1))), max(1,1+round(shift(2))):min(ss(2),ss(2)+round(shift(2))), max(1,1+round(shift(3))):min(ss(3),ss(3)+round(shift(3)))) = false;
        particle(idx) = 0;
    end;
    shifts(i, 1:3) = shift;
    
    

    if (is_wedge_particle)
        particle = real(ifftn(fftn(particle).*wedge));
    end;
    if (sigma_noise > 0)
        if (isempty(mask_for_snr))
            mparticle = mean(particle(:));
            sparticle = std(particle(:), 1);
            particle = (particle-mparticle)/sparticle + sigma_noise * randn(ss);
            particle = (particle * sparticle) + mparticle;
        else
            error('not yet implemented');
        end;
    end;


    
    filename = fullfile(outputdir, ['particle_' strrep(sprintf('%3d', i), ' ', '0') '.em']);
    fprintf(flog, '"%s"    %6.2f %6.2f %6.2f    %5.2f %5.2f %5.2f\n', filename, angle(1), angle(2), angle(3), shift(1), shift(2), shift(3));
    fprintf(fp, '%s\n%s\n', filename, wedge_s);

    tom_emwritec3(filename, tom_emheader(particle));

end;
fclose(fp);


if (exist('append_angles_list', 'var') && numel(append_angles_list)==1 && append_angles_list>0) 
    append_angles_list = init_rotation_list(append_angles_list);
end;
if (exist('append_angles_list', 'var') && size(append_angles_list,2)==3)
    angles = cat(1, angles, append_angles_list);
end;
tom_emwritec3(fullfile(outputdir, 'angles.em'), angles'*pi()/180);
tom_emwritec3(fullfile(outputdir, 'shifts.em'), shifts);

f = fopen(fullfile(outputdir, 'angle_list.txt'), 'w');
fprintf(f, '0 0-%d         %s\n', nparticles-1, fullfile(outputdir, 'angles.em'));
fclose(f);


disp(['FINISHED: ' datestr(now)]);
fprintf(flog, ['\n#finished' datestr(now) '\n']);

fclose(flog);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generates (naively) a random rotation matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angles = init_rotation_list(n)

angles = nan(n, 3);

xlist = nan(3,n);
x = [1,0,0]';
if (1)
    for (i=1:n) 
        R = rand_rot();
        angles(i, :) = tom_rotmatrix2angles(R);
        xlist(:,i) = R * x;
    end;
else
    increment = ceil((360*360*180/n)^(1/3));
    i = 1;
    for (ipsi=1:increment:360)
        for (iphi=1:increment:360)
            for (itheta=1:increment:180)
                angles(i, :) = [ipsi, iphi, itheta];
                [dummy1, dummy2, R] = tom_sum_rotation(angles(i,:), [0,0,0]);
                xlist(:,i) = R * x;
                i = i+1;
            end;
        end;
    end;
end;

%figure(1);
%plot3(xlist(1,:), xlist(2,:), xlist(3,:), '.'); axis equal;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generates (naively) a random rotation matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rot = rand_rot()


n = 5 + ceil(rand()*10);

rot = eye(3);

which_axis = ceil(rand(1,n)*3);
which_axis(1:3) = 1:3;

for (i=1:n)
    alpha = rand()*2*pi();
    switch (which_axis(i))
        case 1
            rot = [[          1,           0,           0]; ...
                   [          0,  cos(alpha), -sin(alpha)]; ...
                   [          0,  sin(alpha),  cos(alpha)]] * rot; 
        case 2
            rot = [[ cos(alpha),           0,  sin(alpha)]; ...
                   [          0,           1,           0]; ...
                   [-sin(alpha),           0,  cos(alpha)]] * rot; 
        case 3
            rot = [[ cos(alpha), -sin(alpha),           0]; ...
                   [ sin(alpha),  cos(alpha),           0]; ...
                   [          0,           0,           1]] * rot; 
    end;
end;






