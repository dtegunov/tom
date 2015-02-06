function ctf=tom_av2_xmipp_plot_ctf_param(f_ctf_param,img_size)
%  TOM_AV2_XMIPP_PLOT_CTF_PARAM plots xmipp ctf-param file
%  
%     tom_av2_xmipp_plot_ctf_param(f_ctf_param)
%  
%    TOM_AV2_XMIPP_PLOT_CTF_PARAM plots xmipp ctf-param file as profile
%                                 useful for comparing ctf-models
%    
%  
%  PARAMETERS
%  
%    INPUT
%     f_ctf_param     filename of ctf-parma file
%     img_size        img_size in pix                     
%
%    OUTPUT 
%     ctf             ctf profile  
% 
%
%  EXAMPLE
%    
%   ctf=tom_av2_xmipp_plot_ctf_param('ctf_4975784_6.ctfparam',4096);
%     
%  
%  REFERENCES
%  
%
%  SEE ALSO
%
%   tom_av2_xmipp_create_ctf_file
%  
%     created by fb ...ole !!
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

disp('...generating ctf profile');

unix(['xmipp_ctf_profile -i ' f_ctf_param ' > '   f_ctf_param '.plot']);

disp('CTF Params:');


data=importdata([f_ctf_param '.plot'],' ',4);
data=data.data;

% transform x-axis 2 pixels

ny_pix=round(img_size./2);
ny_freq=data(end,1);

x_pix=(data(:,2)./ny_freq).*ny_pix;

%figure; plot(x_pix,data(:,4),'b-');

ctf=imresize(data(:,4),[ny_pix 1]);
disp(' ');