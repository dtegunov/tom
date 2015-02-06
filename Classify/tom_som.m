function [result]=tom_som(input_vect,grid_dim,iterations,hist_flag)
%TOM_SOM trains a som for classification and classifies a set of observatiohns using a som
%
%   [result]=tom_som(input_vect,grid_dim,iterations,hist_flag)
%   [classes codevect histogram net]=tom_som(input_vect,grid_dim,iterations)
%
%PARAMETERS
%
%  INPUT
%   input_vect          ...
%   grid_dim            ...
%   iterations          ...
%   hist_flag           ...
%  
%  OUTPUT
%   result      		...
%
%  EXAMPLE
%   im=tom_reshape_stack('/fs/bmsan/apps/tom_dev/data/kerdensom/2d/knispel.em','',2);
%   [classes codevect histogram net]=tom_som(im,[4 4],100);
%   tom_display_som(classes,im,[24 24],codevect);
%
%  REFERENCES
%
%  SEE ALSO
%   tom_kerdendsom tom_display_som
%
%   created by Isar Beder 10/11/06
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

error(nargchk(0, 4, nargin, 'struct'))

sz_input=size(input_vect);
cl_iterations=1;

%generate net
vals=zeros(sz_input(2),2);
vals(:,1)=0;
vals(:,2)=1;

net = newsom(vals,grid_dim);

if (hist_flag==1)
    net.trainFcn='tom_trainr';
    cl_iterations=iterations;
end;


%adjust parameters for training
disp('quiete please ... net is training!');
net.trainParam.epochs = iterations;
net = train(net,input_vect');


for ii=1:cl_iterations

    if (hist_flag==1)
        net.IW{:,:}=net.userdata.IW_hist{:,:,ii};
    end;

    disp('classify particles');
    classes=zeros(sz_input(1),1);
    histogram=zeros(grid_dim(1).*grid_dim(2),1);

    for i=1:sz_input(1)
        tmp=input_vect(i,:)';
        tmp=sim(net,tmp);
        classes(i)=find(str2num((num2str(tmp))));
        histogram(classes(i))=histogram(classes(i))+1;
    end;

    codevect=net.IW{1};
    result(ii).outputMap = codevect';
    result(ii).histogram = histogram;
    result(ii).assignVtoX = classes;
    result(ii).net=net;
end;




