function [outstack,align,histogram] = tom_kerdensom(instructfile, infile, gridsize, binning, maskstruct, normflag)
%TOM_KERDENSOM creates ...
%
%   [outstack,align,histogram] = tom_kerdensom(instructfile, infile,gridsize, binning, maskstruct, normflag)
%
%PARAMETERS
%
%  INPUT
%   instructfile        ...
%   infile              ...
%   gridsize            ...
%   binning             ...
%   maskstruct          ...
%   normflag            ...
%  
%  OUTPUT
%   outstack    		...
%   align       		...
%   histogram    		...
%
%EXAMPLE
%   ... = tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%
% A Novel Neural Network Technique for Analysis and Classification of EM Single-Particle Images
% A. Pascual-Montano, L. E. Donate, M. Valle, M. Bárcena, R. D. Pascual-Marqui, J. M. Carazo 
% Journal of Structural Biology, Vol. 133, No. 2/3, Feb 2001, pp. 233-245 
% http://xmipp.cnb.csic.es/NewXmipp/Applications/Src/KerDenSOM/Help/kerdensom.html%   ...
%
%   created by ... (author date)
%   updated by ...
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


if nargin < 6
    normflag = 1;
end
if nargin < 5
    in_struct.mask.types = {'sphere','rectangle'};
     maskstruct = tom_filtergui('mask',in_struct);
end
if nargin < 4
    binning = 0;
end
if nargin < 3
    gridsize = [5 5];
end

if nargin < 2
    infile = [];
end

outfilename = '/fs/sally/pool-baumeister/tmp/kerdensom/out.dat';
outdir = '/fs/sally/pool-baumeister/tmp/kerdensom';


try;s = load(instructfile);end;
if isequal(instructfile,'2d') || isfield(s,'align2d') 
    %2d data
    header = tom_reademheader(infile);
    dims = header.Header.Size;
    data = tom_reshape_stack(infile,'',binning,'',0,maskstruct);
    create_file(outfilename,data);
    dimflag = 2;
 else
     %3d data
     align = tom_av3_align_sum(s.Align);
     for i=1:size(align,2)
         align(1,i).Shift.X = align(1,i).Shift.X ./ 2^binning;
         align(1,i).Shift.Y = align(1,i).Shift.Y ./ 2^binning;        
         align(1,i).Shift.Z = align(1,i).Shift.Z ./ 2^binning;
     end
    
    filecell = {};
    for i=1:size(align,2)
        filecell{i} = align(i).Filename;
    end
    
    header = tom_reademheader(filecell{1});
    Header = header.Header;
    dims = header.Header.Size;
    in_struct.mask.types = {'sphere3d','cylinder3d'};
    in_struct.rotmask.types = {'sphere3d','cylinder3d'};
    maskstruct = tom_filtergui('mask',in_struct);
    
    data = tom_reshape_stack(filecell,align,binning,'',0,maskstruct);
    
    %create_file(outfilename,data);    
    dimflag = 3;
end

if normflag == 1
%    data = normalize(data);
end

inputVectors = data';
outputMapHeight = gridsize(1);
outputMapWidth = gridsize(2);
initialReg = 1000;
finalReg = 100;
regSteps = 10;
som_Nsteps = 200;
initialMap = [];

[outputMap, alpha, U, histogram, assignVtoX, assignXtoV] = kerdensom(inputVectors, outputMapHeight,  outputMapWidth, initialReg,  finalReg, regSteps, som_Nsteps, initialMap,'norm');
%outstack = zeros(dims(1),dims(2),gridsize(1).*gridsize(2));
if dimflag == 2
    outstack=reshape(outputMap,[dims(1)./ 2^binning dims(2)./ 2^binning gridsize(1).*gridsize(2)]);
else
    outstack=reshape(outputMap,[dims(1)./ 2^binning dims(2)./ 2^binning dims(3)./ 2^binning gridsize(1).*gridsize(2)]);
end

if dimflag == 3
    
    for i=1:size(assignVtoX,1)
        align(i,end).Class = assignVtoX(i);
    end

else
    %todo 2d alignment file
end

%unix(['cd ' outdir ' ; setenv LD_LIBRARY_PATH /usr/lib ; /usr/local/apps/Xmipp-0.9/bin/xmipp_kerdensom -din ' outfilename ' -cout test_kerdensom -xdim ' num2str(gridsize(1)) ' -ydim ' num2str(gridsize(2)) ' -verb 1 -norm']);


%process results
%outstack = codfile2stack(outdir,dims,gridsize);
%histogram = read_histogramfile(outdir,gridsize);
%align='';

function stack = normalize(stack)

mea = mean2(stack);
sd = std2(stack);

stack = (stack-mea)./sd;


function image = codfile2stack(directory,dims,gridsize)

fid = fopen([directory '/test_kerdensom.cod'],'rt');

%read file header
tline = fgets(fid);

image = zeros(dims(1),dims(2),gridsize(1).*gridsize(2));

for i=1:gridsize(1).*gridsize(2)
    tline = fgets(fid);
    tmp =  sscanf(tline,'%f',[dims(1) dims(2)]);
    image(:,:,i) = tmp;
end

fclose(fid);

function histogram = read_histogramfile(directory,gridsize)

fid = fopen([directory '/test_kerdensom.his'],'rt');

%read file header
tline = fgets(fid);

histogram = zeros(gridsize(1).*gridsize(2),1);

for i=1:gridsize(1).*gridsize(2)
    tline = fgets(fid);
    tmp =  sscanf(tline,'%f',2);
    histogram(i) = tmp(2);
end

fclose(fid);



function image = codfile2stack3d(directory,dims,gridsize)

fid = fopen([directory '/test_kerdensom.cod'],'rt');

%read file header
tline = fgets(fid);

image = zeros(dims(1),dims(2),dims(3),gridsize(1).*gridsize(2));

for i=1:gridsize(1).*gridsize(2)
    tline = fgets(fid);
    tmp =  sscanf(tline,'%f',[dims(1) dims(2) dims(3)]);
    image(:,:,:,i) = tmp;
end

fclose(fid);


function create_file(filename,data)

fid = fopen(filename,'wt');

fprintf(fid,'%i %i\n',size(data,2),size(data,1));

for i=1:size(data,1)
    for j=1:size(data,2)
        fprintf(fid,'%0.5f ', data(i,j));
    end
    fprintf(fid,' particle_%i\n',i);
end

fclose(fid);


