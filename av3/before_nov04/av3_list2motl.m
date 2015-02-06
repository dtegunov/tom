function motl = av3_list2motl(list)
%
%   motl = av3_list2motl(list)
% 
%   converts strange format list from tom_particles to motl. For MOTL format
%   see tom_picker or tom_chooser. NOTE: The format to proceed processing
%   data with the AV3 package is the motl format. 
%
% FF 11/03/03

nparts = size(list,2);
motl = zeros(20,nparts);
motl(4,:) = list(1,:);
motl(8:10,:) = list(6:8,:);

