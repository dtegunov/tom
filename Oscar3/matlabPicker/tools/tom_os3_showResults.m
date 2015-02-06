function res = my_showResults(peaks, img, template, threshold)

if(my_fftSelect(img) == 3)
    res = img;
    return;
end;


mx = my_max(peaks);

bool = peaks >= mx*threshold;
bool2 =zeros(size(img));
mask = ones(size(template));
for x = 1:size(peaks,1)
    for y=1:size(peaks,2)

        if(bool(x,y))
            
            bool2 = tom_paste(bool2,mask,[x-floor(size(mask,1)/2)+1 y-floor(size(mask,2)/2)+1]);
            
        end;
    end;
end;

%my_show(template);

erod = imerode(bool2,[1 1 1; 1 1 1; 1 1 1]);


res = img + mean(mean(img))*(bool2-erod);
%imtool(res)

