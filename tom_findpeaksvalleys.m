function [ peaks, valleys ] = findpeaksvalleys( data, searchdistance )

peaks = [];
valleys = [];

if size(data,1)<size(data,2)
    data = data';
end;

for x=1:size(data,1)
    segment = data(max(1, x - searchdistance):min(size(data, 1), x + searchdistance));
    if data(x) == max(segment)
        peaks = [peaks; x, data(x)];
    elseif data(x) == min(segment)
        valleys = [valleys; x, data(x)];
    end;
end;

end

