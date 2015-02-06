%%
%determine the maximum of the peaks
function [v p] = tom_os3_peakMax(peaks,psr)

%%
%get peak maximum
%[peakMax.val peakMax.pos]  = tom_os3_max(peaks);
[v p]  = tom_os3_max(peaks);
return;

%get psr maximum
[psrMax.val psrMax.pos] = tom_os3_max(psr);

%if the distance between peak maximum and psr maximum is close, use this
%peak
if(sqrt((peakMax.pos - psrMax.pos) * (peakMax.pos - psrMax.pos)') <= 10)
    
    v = peakMax.val; 
    p = peakMax.pos;

else        
%%
%else use linear combination of peaks and psr as criterion
    if(peakMax.val > 0)
        [v p] = tom_os3_max(peaks .* psr);
    else
        v = peakMax.val; 
        p = peakMax.pos;
    end;

end;
%%
%set return values
v = peakMax.val; 
p = peakMax.pos;