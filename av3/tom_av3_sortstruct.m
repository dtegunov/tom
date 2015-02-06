function structure = tom_av3_sortstruct(structure, fieldname, order, step)



if step == 0
    step = 1;
end

if nargin < 3
    order = 'ascend';
end

sortmatrix = zeros(1,size(structure,2));

for i=1:size(structure,2)
    sortmatrix(i) = structure(step,i).(fieldname);
end

[sorted,idx] = sort(sortmatrix,2,order);

for i=1:size(structure,2)
    for j=1:size(structure,1)
        outstruct(j,i) = structure(j,idx(i)); 
    end
end

structure = outstruct;