function docf=tom_av2_xmipp_bin_doc(docf,fact,output_doc,find_what,replace_with)
%TOM_AV2_XMIPP_BIN_DOC bins doc file (shifts are multiplied by a factor)
%
%   doc_out=tom_av2_xmipp_bin_doc(input_doc,fact,output_doc)
%
%  
%
%PARAMETERS
%
%  INPUT
%   docf            filename of the input doc
%   fact            factor shifts should be multiplied (bin1==0.5) (bin2==0.25) 
%   output_doc      filename of output doc
%   find_what       string 2 be replaced 
%   replace_with    replacement string
%
%  OUTPUT
%   doc_out         binned doc in memory
%
%EXAMPLE
%  
%   doc_out=tom_av2_xmipp_bin_doc('holo_r7i14_chi64.doc',4,'holo_r7i14_chi256.doc','/parts_holo64/','/parts_holo256/');
%   
%
%REFERENCES
%
%SEE ALSO
%   
%   by fb 
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


if (ischar(docf))
    docf=tom_xmippdocread(docf);
end;



for i=1:length(docf)
    if (isempty(find_what)==0)
        tmp=strrep(docf(i).name,find_what,replace_with);
        if (strcmp(tmp,docf(i).name)==1)
            disp(['warning repl. does not work: ']);
            disp([tmp ' == ' docf(i).filename]);
        else
            docf(i).name=tmp;
        end;
     end;
    
    docf(i).xoff=docf(i).xoff.*fact;
    docf(i).yoff=docf(i).yoff.*fact;
end;

if (isempty(output_doc)==0)
    tom_xmippdocwrite(output_doc,docf);
end;

