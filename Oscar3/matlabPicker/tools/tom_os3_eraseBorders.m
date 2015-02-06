function res = tom_os3_eraseBorders(volume,template,binning)


    res = zeros(size(volume),'single');

    innerVolume = ones(size(volume)-(size(template)+1),'single');
    
    res = tom_os3_pasteCenter(res,innerVolume);
    
    res = volume .* res;
    