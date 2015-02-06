function tom_createctfstack(directory,size_out,file_out,binning)

if nargin < 4
    binning = 1;
end

[dircell, flagcell_out] = get_dircontents(directory, {}, {});

%create new stack
tom_emwritec(file_out,[size_out,size_out,length(dircell)],'new');

%create aligment file for metadata
align2d = tom_av2_create_alignfromstack(file_out);


h = waitbar(0,'Getting directory list...');
for i=1:length(dircell)
    file = dircell{i};
    image = tom_emreadc(file,'binning',binning);
    image.Value =single(image.Value);
    image.Value = tom_norm(log(abs(tom_ps(image.Value))),'3std');
    ctf_out = tom_create_ctf(image.Header.FocusIncrement./10000, image.Value, image.Header.Objectpixelsize/10, image.Header.Voltage/1000, image.Header.Cs);
    image.Value(1:size(image.Value,1)./2,:) = tom_norm(ctf_out(1:size(image.Value,1)./2,:).^2,'2std');
    image.Value = imresize(image.Value,[size_out,size_out],'bicubic');
    tom_emwritec(file_out,image.Value,'subregion',[1 1 i],[size_out,size_out,1]);
    align2d(1,i).filename = file;
    align2d(1,i).ccc = image.Header.FocusIncrement./10000;
    align2d(1,i).class = image.Header.Defocus./10000;
    waitbar(i./length(dircell),h,[num2str(i), ' of ', num2str(length(dircell)), ' files processed.']);
end

%save the metadata file
[pathstr, name] = fileparts(file_out);
save([pathstr '/' name '.mat'],'align2d');

 close(h);   
    

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
