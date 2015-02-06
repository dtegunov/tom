function [st,Align]=tom_pytomXmlPartlist_read(filename,atribute)
%tom_pytomXmlPartlist_read reads xml from pytom in a matlab struct
%
%   st=tom_pytomXmlPartlist_read(filename)
%
%PARAMETERS
%
%  INPUT
%   filename            input filename
%   atribute            ('Particle') attribute 4 parsing   
%  
%  OUTPUT
%   st                 output struct
%   Align              ouput struct      
%
%EXAMPLE
%   read image
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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


%parse inputs

if (nargin < 2)
    atribute='Particle';
end;

tmpStruct = parseXML(filename);
tmpStruct=tmpStruct.Children(2).Children(2);

%parse number of particles
zz=0;
for i=1:length(tmpStruct.Children)
    if (strcmp(tmpStruct.Children(i).Name,atribute))
        zz=zz+1;
        idx(zz)=i;
    end;
end;
disp(['Found ' num2str(zz) ' entres found']);




for i=1:length(idx)
    st(i).Name=tmpStruct.Children(idx(i)).Attributes.Value;
    for ii=1:length(tmpStruct.Children(idx(i)).Children)
        if (strcmp(tmpStruct.Children(idx(i)).Children(ii).Name,'#text')==0)
            zz2=0;
            for iii=1:length(tmpStruct.Children(idx(i)).Children(ii).Attributes)
                if (isnan(str2double(tmpStruct.Children(idx(i)).Children(ii).Attributes(iii).Value))==0)
                    zz2=zz2+1;
                    tmp_vect(zz2)=str2double(tmpStruct.Children(idx(i)).Children(ii).Attributes(iii).Value);
                end;
            end;
            
            if (strcmp(tmpStruct.Children(idx(i)).Children(ii).Name,'Rotation'))
                new=[tmp_vect(2) tmp_vect(1) tmp_vect(3)];
                tmp_vect=new;
            end;
            call=['st(i).' tmpStruct.Children(idx(i)).Children(ii).Name '=tmp_vect;'];
            eval(call);
            clear('tmp_vect');
        end;
    end;
    if (mod(i,100)==0)
        disp([num2str(i) ' entries processed!']);
    end;
end;


for i=1:length(st)
    Align(i).Filename=st(i).Name;
    Align(i).Shift.X = st(i).Shift(1);
    Align(i).Shift.Y = st(i).Shift(2);
    Align(i).Shift.Z = st(i).Shift(3);
    Align(i).Angle.Phi = st(i).Rotation(1);
    Align(i).Angle.Psi = st(i).Rotation(3);
    Align(i).Angle.Theta = st(i).Rotation(2);
    Align(i).Angle.Rotmatrix = [];
    Align(i).CCC=st(i).Score(1);
    Align(i).Class=0;
    Align(i).Tomogram.AngleMin=-60;
    Align(i).Tomogram.AngleMax=60;
    Align(i).Tomogram.Header.Size=[64 64 64];
    Align(i).Tomogram.Offset = 0;
    Align(i).Tomogram.Binning = 0;
    Align(i).Tomogram.Regfile = '';
    Align(i).Tomogram.Position.X =st(i).PickPosition(1); %Position of particle in Tomogram (values are unbinned)
    Align(i).Tomogram.Position.Y =st(i).PickPosition(2);
    Align(i).Tomogram.Position.Z =st(i).PickPosition(3);
    Align(i).ProjectionClass = 0;
    Align(i).NormFlag = 0; %is particle phase normalized?
    Align(i).Filter = [0 0]; %is particle filtered with bandpass?

end;



function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems
% with very deeply nested trees.
try
    theStruct = parseChildNodes(tree);
catch
    error('Unable to parse XML file %s.',filename);
end


% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1, numChildNodes);
    
    children = struct(             ...
        'Name', allocCell, 'Attributes', allocCell,    ...
        'Data', allocCell, 'Children', allocCell);
    
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
    'Name', char(theNode.getNodeName),       ...
    'Attributes', parseAttributes(theNode),  ...
    'Data', '',                              ...
    'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end









