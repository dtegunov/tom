function classes = tom_av2_spiderclasses(refname,stackname,paramsfile)

refs = tom_readspiderstack(refname);
stack = tom_readspiderstack(stackname);

A = importdata(paramsfile);

%allocate classes stack
classes = zeros(stack.Header.Size(1), stack.Header.Size(2), refs.Header.Size(3).*2,'single');

%fill in references
j = 1;
for i=1:refs.Header.Size(3)
    classes(:,:,j) = tom_norm(refs.Value(:,:,i),'mean0+1std');
    j = j + 2;
end

%create class averages
for i=1:max(A.data(:,6))
    indx = find(A.data(:,6) == i);

    if ~isempty(indx)
        part = zeros(stack.Header.Size(1), stack.Header.Size(2),'single');
        for j=1:size(indx,1)
            if A.data(indx(j),17) == 1
                part = part + stack.Value(:,:,indx(j));
            else
                part = part + tom_mirror(stack.Value(:,:,indx(j)),'x');
            end
        end

        part = part ./ size(indx,2);
        classes(:,:,i.*2) = tom_norm(part,'mean0+1std');
    end
    
end

figure;tom_dspcub(classes);