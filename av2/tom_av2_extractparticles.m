function align=tom_av2_extractparticles(ccfmap, anglesmap, templatemap, tsize, number, filename, binning, align)
%TOM_AV2_EXTRACTPARTICLES creates ...
%
%   align=tom_av2_extractparticles(ccfmap, anglesmap, templatemap, tsize, number, filename, binning, align)
%
%PARAMETERS
%
%  INPUT
%   ccfmap              ...
%   anglesmap           ...
%   templatemap         ...
%   tsize               ...
%   number              ...
%   filename            ...
%   binning             ...
%   align               ...
%  
%  OUTPUT
%   align               ...
%
%EXAMPLE
%   ... = tom_av2_extractparticles(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

binning = 2^binning;

if nargin == 6
    align = struct();
    start = 1;
else
    if size(align,2) ~= 1
        start = size(align,2)+1;
    else
        start = 1;
    end
end


%load trainingdata
%t = load('/home/korinek/trainingdata.mat');
%t = load('/fs/sally/pool-baumeister/data/170305/log/pick/trainingdata.mat');
%svmStruct = svmtrain(double(t.trainingdata),t.g');
%lauf2 = 1;
%lauf = start;
for lauf=start:start+number-1
%while lauf < start+number
[pos val ccfmap] = tom_peakc(ccfmap,tsize./2);
 %   radius = round(tsize./2);
 %   im = tom_emreadc(filename,'subregion',[round(pos(1))-2*radius round(pos(2))-2*radius 1],[4*radius-1 4*radius-1 0]);
  %  im = reshape(double(im.Value),[],1);
    %class = classify(tom_norm(im',1),tom_norm(double(t.trainingdata),1),t.g,'diaglinear');
    
   % Group = svmclassify(svmStruct, im');
   % Group
    
    %if Group == 1
        align(1,lauf).dataset = '';
        align(1,lauf).filename = filename;
        align(1,lauf).position.x = round(pos(1).*binning);
        align(1,lauf).position.y = round(pos(2).*binning);
        align(1,lauf).class = ['auto_picked_' num2str(templatemap(round(pos(1)),round(pos(2))))];
        align(1,lauf).ref_class = templatemap(round(pos(1)),round(pos(2)));
        align(1,lauf).radius = tsize.*binning./2;
        align(1,lauf).color = [1 0 0];
        align(1,lauf).shift.x = 0;
        align(1,lauf).shift.y = 0;
        align(1,lauf).angle = -anglesmap(round(pos(1)),round(pos(2)));
        align(1,lauf).isaligned = 1;
        align(1,lauf).ccc = val;
        align(1,lauf).quality = 0;
        align(1,lauf).normed = 'none';
        lauf=lauf+1;
    %end
    %lauf2 = lauf2+1;
    %if lauf2 > 300
    %    disp('>300');
    %    return;
    %end
end
%disp('image done');
