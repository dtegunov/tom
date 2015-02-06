function tom_os3_saveResults(peakList, res, resultNo,flags)

vol = res.img{resultNo};
template = res.templates{resultNo};

dim = tom_os3_fftSelect(vol);


%%
%choose folder where to save the picklist images
folderName= inputdlg('Select name of new folder','Name of new folder');
if(length(folderName) == 0)
    folderName = '';
end;
folderName = ['/' char(folderName)];

directoryName = uigetdir('.', 'Pick a directory');
mkdir(directoryName,char(folderName));

counter = 1;

%%
%for each picklist entry do:
for i = 1:length(peakList)
    
    peak = peakList{i};
%positions of peak (center)
    x = peak(1);
    y = peak(2);
    z = peak(3);
%calculated value    
    value = peak(4);
%angle of best match    
    angle = peak(5);        

%get the subvolume at pos x,y,z with the size of the template    
    subvolume = tom_os3_subVolume([x y z],vol,template,'center');
    
%save the subvolume    
    subvolume= tom_emheader(subvolume);
    tom_emwrite([directoryName folderName num2str(counter) '_pick' '_angle_' mat2str(angle) '.em'],subvolume);
    
    counter = counter +1;
end;

%%
%save picklist itself

