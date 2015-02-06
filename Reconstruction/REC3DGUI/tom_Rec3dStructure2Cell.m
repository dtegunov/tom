function CellOut = tom_Rec3dStructure2Cell(StructureIn)
%TOM_REC3DSTRUCTURE2CELL is a module of TOM_REC3DGUI.
%
%   CellOut = tom_Rec3dStructure2Cell(StructureIn)
%
%   It converts a structure to a cell array which can be displayed 
%   in a listbox.
%
%PARAMETERS
%
%  INPUT
%   StructureIn        MATLAB structure
%  
%  OUTPUT
%   CellOut            MATLAB cell array
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%
%   01/01/07 ME
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


% Write fieldnames to cell
CellOut = fieldnames(StructureIn);

% Add fieldvalues
for k = 1:size(CellOut,1)
    
    FieldValue = getfield(StructureIn, CellOut{k});
    
    if isstruct(FieldValue) == 1
         CellOut{k} = [CellOut{k} '=[<1x1 struct>]'];
    else
         CellOut{k} = [CellOut{k} '=[' num2str(FieldValue) ']'];
    end
    
end

