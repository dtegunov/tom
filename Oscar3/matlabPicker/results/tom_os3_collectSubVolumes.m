%%
function subVolumes = tom_os3_collectSubVolumes(picklist,volume)

%%  if volume is numeric, continue, else load the specified file
    if(ischar(volume))
    
        volume = tom_emread(volume);
        volume = single(volume.Value);
        
    end;
    
%%  pick each element of the picklist out of the volume
    subVolumes = {};
    for i = 1:length(picklist)
        
        pick = picklist{i};
        
        object = tom_os3_subVolume(pick(1:3),volume,zeros(pick(6:end),'single'),'center');
        
        subVolume.volume = object;
        subVolume.pick   = pick;
        subVolumes{length(subVolumes)+1} = subVolume;

    end;