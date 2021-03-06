function [res options]= tom_os3_modifyImage(image,options)
%tom_os3_modifyImage
%
%   tom_os3_modifyImage(image,options)
%   applies bandpassFiltering to the image
%   applies a point spread function to the image
%   
%   the modifications are specified in the options structure under the
%   attribute  modifications
%PARAMETERS
%
%  INPUT
%   image    - the image to be modified
%   options  - option structure generated by tom_os3_readOptions
%  
%  OUTPUT
%   res      - the resulting image
%   options  - the options structure. usefull if it has been modified
%              during the execution of this function
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_os3_readOptions
%
%   created by TH2 07/07/07
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

%%  apply bandpass filtering if lower frequency value reasonable
    if(options.modifications.bandpass.low > 0)
        
        lower   = options.midifications.bandpass.low;
        upper   = options.midifications.bandpass.high;
        
        image = tom_bandpass(image,lowerBand,upperBand,0);
        
    end;
%%  apply point spread function if defined
    if(ischar(options.psf.file) && ~strcmp(options.psf.file,'none'))
    %   load psf file to mem
        try
            psf                  = tom_emread(options.psf.file);
        catch
            err     = lasterror;
            error(['The specified psf coused an error. Please check your option file again!' err.message]);
        end;
        
        options.psf.file     = psf.Value;
        
    end;    
    
    if(isnumeric(options.psf.file) && ~ischar(options.psf.file))    
        
        if(isequal(size(image),size(options.psf.file)) )
            image   = tom_apply_weight_function(image,options.psf.file);
        else
            img     = tom_os3_pasteCenter(zeros(size(options.psf.file)),image);
            img     = tom_apply_weight_function(img,options.psf.file);
            image   = tom_os3_subVolume(size(options.psf.file)/2+1,img,image,'center');
            
        end;
    end;
    
    res     = image;