function tom_av2_xmipp_filter_som(align2d,output_doc,output_sel)
%  tom_av2_xmipp_filter_som generates a filtered doc, and sel according to
%  the selection made in tom_av2_stackbrowser
%  
%     tom_av2_xmipp_filter_som(align2d,output_doc,output_sel)
%  
%  PARAMETERS
%  
%    INPUT
%     align2d        filtered align2d form tom_av2_xmipp_somvs2stack and
%                    tom_av2_stackbrowser
%     output_doc     doc file of the selected particles
%     output_sel     sel file of the selected particles 
%     
%
%  EXAMPLE
%  
%  tom_av2_xmipp_filter_som('filt.mat','filt.doc','filt.sel');
%
%  for re-using in xmipp kerdensom 
%  mkdir ML2d/dum_ml1
%  cp SOM/ML2ref_dragon_mbig5/out_som/filt2.sel ML2D/dum_ml1/ml2d_ref000004.sel
%  cp SOM/ML2ref_dragon_mbig5/out_som/filt2.doc ML2D/dum_ml1/ml2d_it000001.doc
% 
%  %enter dum_ml1 in "Directroy where you have previously run ML2d" of xmipp kerdensom gui  
%  %enter 4  in "The number of the class to use"  
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_av2_xmipp_somvs2stack,tom_av2_stackbrowser
%  
%     created by FB 03/02/10
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if isstruct(align2d)==0
    load(align2d);
end;
   

[a b]=unix(['cat ' align2d(1,1).filename  '.doc > ' output_doc '.tmp']);
if (a==1)
    error(b);
end;
[a b]=unix(['cat ' align2d(1,1).filename  '.sel > ' output_sel '.tmp']);
if (a==1)
    error(b);
end;

for i=2:size(align2d,2)
    unix(['cp ' align2d(1,i).filename  '.doc in.txt.tmp']);
    unix(['sed -i ''1d''' ' in.txt.tmp']); 
    unix(['cat ' output_doc '.tmp in.txt.tmp > ' output_doc]);
    
    unix(['cp ' output_doc ' ' output_doc '.tmp']);
    
    unix(['cat ' output_sel '.tmp ' align2d(1,i).filename  '.sel > ' output_sel]);
    unix(['cp ' output_sel ' ' output_sel '.tmp']);
    disp(['processing class nr: ' num2str(i) ' of ' num2str(size(align2d,2)) ]);
end;

unix(['rm ' output_sel '.tmp']);
unix(['rm ' output_doc '.tmp']);
unix(['rm in.txt.tmp']);




