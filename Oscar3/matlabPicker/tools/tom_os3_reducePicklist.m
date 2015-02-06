function [newPickList newStack] = tom_os3_reducePicklist(picklist,goodFlags,particleStack)

stackExists = exist('particleStack');

newPickList = {};
newStack = [];
for i=1:length(goodFlags)
    
    if(goodFlags(i) == 1)
        newPickList{length(newPickList)+1} = picklist{i};
        if(stackExists)
            newStack(:,:,size(newStack,3)+1) = particleStack(:,:,i);
        end;
    end;
    
end;