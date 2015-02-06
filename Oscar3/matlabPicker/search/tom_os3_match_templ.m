function comb_cc=tom_os3_match_templ(image,templ_dir,options)
%  TOM_OS3_MATCH_TEMPL calculates ccf for given templ directory
%  
%     [pos comb_cc stat]= tom_os3_match_templ(image,templ_dir,options)
%  
%  PARAMETERS
%  
%    INPUT
%     image               search image filename or matrix
%     templ_dir           image to be aligned
%     options             options struct
%
%    
%    OUTPUT
%    
%     comb_cc           final peaks
%
%
%  EXAMPLE
%  
%  options.f_xcf_opt=0.2;
%  options.f_xcf_opt=0;
%  options.f_xcf_opt=0.8;
%  
%  cc=tom_os3_match_templ(tom_bin(im.Value,2),'output/tmpl_test/',options);
%  
%  % to be done!
%  %options.psf='none'; % or filename %or 'fit_st'
%  %options.filter.apply=0; % to switch off
%
%  
%
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by FB 
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


if (nargin<3)
    f_xcf=0.2;
    f_soc=0.8;
    f_psr=0;
else
    f_xcf=options.f_xcf_opt;
    f_soc=options.f_soc_opt;
    f_psr=options.f_psr_opt;
end;


%%  set system wisdom file for fftw
    hostname = tom_os3_getHostname;
   % disp(['Processing ' job.volumeFile ' on: ' hostname]); 
    
%    
%    log=loggerCl;
%    log.filename=[options.result(1).IO.outputdir '/logs/log.txt'];
%    log.entry('Progress:Filename',tom_os3_filename_and_host(job.volumeFile),hostname);
%    
   
   %options.result(1).log.entry('Progress:Filename',tom_os3_filename_and_host(job.volumeFile),hostname);
   




%%  init function


dimension       = size(image);
tmp  = dir([templ_dir '/*.em']);

for i=1:length(tmp) 
    templateListe{i}=[templ_dir '/' tmp(i).name];
end;
   
    



%%  init return values (allocate memory
   
    ccc         = -ones(size(image)); 
    psr         = ccc;
    autocorr    = ccc;
    angles      = ccc;
    
%%  if volume is empty, nothing can be found -> return zeros    
    if(std(image(:)) == 0)
        ccc     = ccc *0;
        psr     = ccc *0;
        autocorr= ccc *0;
        return;
    end;
        




%% fro each template do
for templateIterator =1:length(templateListe)
    
    
    
    template = tom_emreadc3(templateListe{templateIterator});
    template = single(template.Value);
    
    templateSize = size(template);
    if(length(templateSize) < 3)
        templateSize(3) = 0;
    end;
    
    [mask maskSize]= tom_os3_sphereMask(template);
    
    %normalize volume here, -> less calculations in for loop in
    imageMean = tom_os3_mean(image,mask,maskSize);
    imageSTD  = tom_os3_std(image,imageMean,mask,maskSize);
    fimage    = fftn(image);
    calculationAvailable = true;
    
    
    
    
    template = tom_norm(100000*(template+1),'mean0+1std');
    
    [ tmp templateMean templateSTD] = tom_os3_normUnderMask(template,mask);
    
    t   = template.* mask;
    % [t options]  = tom_os3_modifyImage(t,options);
    
    [cccMap options]    = tom_os3_corr(image,t);
    cccMap = single(cccMap);
    
    
    if (f_psr~=0)
        opt.analysis.type   = 'PSR';
        psrMap=single(tom_os3_analyzePeak(image,cccMap,t,opt));
    else
        psrMap=ones(size(cccMap),'single');
    end;
    
    if (f_soc~=0)
        opt.analysis.type   = 'Autocorrelation';
        autocorrMap = single(tom_os3_analyzePeak(0,cccMap,t,opt));
    else
        autocorrMap = ones(size(cccMap),'single');
    end;
    
    
    
   
    opt.analysis.ccc=f_xcf~=0;
    opt.analysis.psr=f_psr~=0;
    opt.analysis.autocorr=f_soc~=0;
    
    
    [ ccc angles psr autocorr] = tom_os3_bestPeak(ccc,cccMap, ...
        angles,ones(size(angles),'single'), ...
        psr,psrMap, ...
        autocorr,autocorrMap, ...
        dimension,opt);
   
    
end;

clear('angles');


%%  set the results
%   this job is not returned
ccc       = single(ccc);
psr       = single(psr);
autoc     = single(autocorr);


comb_cc = f_xcf * ccc + f_psr * psr + f_soc * autoc;



%log.entry('Progress:TemplateMatching',tom_os3_filename_and_host(job.volumeFile),hostname,length(pickList),mean(tmp_cc));

