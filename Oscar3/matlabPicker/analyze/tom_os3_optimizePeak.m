
%%
%determine a better peak value as a function of the anglespace
%set staring and ending angle of the correlation interval and reduce the
%angleIncrement 
%determine the maximum of this more accurate scan
function res = tom_os3_optimizePeak(img,template,flags,currentAngle,mode)

dim = tom_os3_fftSelect(img);

angle = flags.angle;

angleIncrement = angle.angleIncrement;
if(dim == 2 )

    if(strcmp(mode,'peak'))
        
        angle.start = currentAngle - floor(angleIncrement/2);
        angle.end = currentAngle + floor(angleIncrement/2);
        angle.angleIncrement = 1;
        
    elseif(strcmp(mode,'interval'))
        
        angle.start = currentAngle - angleIncrement +1;
        angle.end = currentAngle + angleIncrement - 1;
        angle.angleIncrement = floor(angleIncrement / 5); %2*angInc / 10

    end;
    
else(dim == 3)

end;

flags.angle = angle;

%%
%determine more accurate peak values
% [peakMax angles psr] = tom_os3_findTemplate(img,template,flags);
s.volume      = img;
s.template    = template;
s.angleList   = tom_os3_angleList(flags,dim);
s.angleOffset = 0;
s.coordinates = 0;
s.flags       = flags;  
s.dimension   = dim;
[peakMax angles psr] = tom_os3_findTemplate(s);


%create mask around the peak
if(dim == 2)
    peakMask = ones(10,10,1);
else
    peakMask = ones(10,10,10);
end;

%better- faster approach - cut out the values around the center
peakMask = tom_os3_pasteCenter(zeros(size(img)),peakMask);

peakMax = peakMax .* peakMask;
angles = angles .* peakMask;
psr = psr .* peakMask;


[val pos] = tom_os3_max(peakMax .* psr);

res.val = val;

if(dim == 2)
    res.angle = angles(pos(1),pos(2));
else
    %dim 3 
end;
    