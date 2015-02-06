%%
%returns picklist of the best n(formally known as threshold) peaks
function pickList = tom_os3_returnPicklist(peaks,angles,templateSize,n,options)

%get dimension of the data
dim = options.dimension;

%set the radius of the peak area - this area will be set to zero after one
%of the peaks has been selected
radius = floor(templateSize(1)/2);
pickList = {};

counter = 1;
while(counter <= n)
    
    [coordinates value peaks] = tom_peak(peaks,radius);

    if(size(template,3) <= 1)
        %if is 2d template
        coordinates(3) = 1;
    end;

    if(iscell(angles))
        a = angles.value;
        angl = a(coordinates(1),coordinates(2),coordinates(3));
    else
        angl = angles(coordinates(1),coordinates(2),coordinates(3));
    end;
    
    pickList{counter} = [coordinates value angl templateSize];
    counter = counter + 1;
    
end;