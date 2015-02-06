function tom_ctfcorrect_dir(directory, outdir, method, filefilter, mtf, cutoff, defocus,parallelstruct)

if nargin < 8
    parallelstruct = [];
end

if nargin < 4 || isempty(filefilter)
    [dircell, flagcell] = get_dircontents(directory, {}, {});
else
    s = load(filefilter);
    dircell = s.particlepicker.filelist;
    flagcell = s.particlepicker.filefilter;
end


%parallel execution, TODO
if ~isempty(parallelstruct) && ~isempty(fieldnames(parallelstruct))
    numtasks = parallelstruct.number_of_tasks;
    jm = findResource('jobmanager','name',parallelstruct.jobmanager);
    j = createJob(jm,'Name','autopicker');
    pdd={'/fs/pool/pool-bmsan-apps/tom_dev/Sptrans' '/fs/pool/pool-bmsan-apps/tom_dev/Analysis/' '/fs/pool/pool-bmsan-apps/tom_dev/IOfun/' '/fs/pool/pool-bmsan-apps/tom_dev/av2/' '/fs/pool/pool-bmsan-apps/tom_dev/Filtrans/' '/fs/pool/pool-bmsan-apps/tom_dev/Geom/'};
    set(j,'FileDependencies',pdd);
    packages=tom_calc_packages(numtasks,size(imagecell,2));
    
    disp('Sending pick jobs to workers...');
    numtasks2 = numtasks;
    for i=1:numtasks
    end
    
%local execution
else
    for i=1:length(dircell)
        if flagcell{i} == 1 && tom_isemfile([outdir '/' dircell{i}]) ~= 1;
            file = dircell{i};
            img = tom_emreadc([directory '/' file]);

            if nargin>6
                %img.Header.FocusIncrement = defocus.*10000;
            end
            if nargin > 4
                img = tom_mtfdeconv(img,method,mtf,cutoff);
            else
                img = tom_mtfdeconv(img,method);
            end
            img.Header.EM.Magic(4) = 5;
            img.Value = single(img.Value);
            tom_emwritec([outdir '/' file],img);

            disp(file);
        end
        

    end

    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell, flagcell_out] = get_dircontents(directory, dircell, flagcell)

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Getting directory list...');

dirlist = dir(directory);
files = size(dirlist,1);
j = size(dircell,2);
if j == 0
    j = 1;
end

%sort out all the em-files
for i = 1:files
    if dirlist(i).isdir == 0
        if isempty(strmatch(dirlist(i).name,dircell))
            if tom_isemfile([directory '/' dirlist(i).name]) == 1
                dircell{j} = dirlist(i).name;
                try
                    flagcell_out{j} = flagcell{j};
                catch
                    flagcell_out{j} = 1;
                end
                j = size(dircell,2) + 1;
            end
        end
    end
    waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end

set(findobj('Tag','fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);
