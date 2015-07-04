function [ output ] = tom_mavg( data, window )

windowhalf = window / 2;
if size(data, 2) > size(data, 1)
    data = data';
end;
output = data;

for i=1:size(data, 1);
    output(i) = mean(data(max(1, i-windowhalf):min(size(data,1), i + windowhalf)));
end