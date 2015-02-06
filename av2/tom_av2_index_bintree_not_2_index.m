function index=tom_av2_index_bintree_not_2_index(bin_tree_not,base)
%TOM_AV2_INDEX_BINTREE_NOT_2_INDEX calculates stack positon out of left
%right notation ...0001100
%
%    index=tom_av2_index_bintree_not_2_index(bin_tree_not,base)
%PARAMETERS
%
%  INPUT
%   bin_tree_not     b-tree notation bayer-man rocks 
%   base             number of classes 
%   
%  OUTPUT
%   index    unique stack position                 
%
%EXAMPLE
%     index=tom_av2_index_bintree_not_2_index('100001101',2) % for a binary tree
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_index_search, 
%   tom_av2_index_plot_indexstack.m
%
%   created by fb (eckster)
%   updated by fb (eckster) adapeted for any base needed for n-classes!
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



ideg=base2dec(bin_tree_not,base);
num_tmp=0;
for ii=1:size(bin_tree_not,2)-1
    num_tmp=num_tmp+base^(ii);
end;
index=num_tmp+ideg+1;