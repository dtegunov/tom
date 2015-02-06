function im_out=tom_apply_filter(im_in,filter_st)
%TOM_APPLY_FILTER performs image filtering according to the given structure
%
%   im_out=tom_apply_filter(im_in,filter_st)
%
%PARAMETERS
%
%  INPUT
%   im_in               reference image
%   filter_st           filter structure:
%                                  filter.Apply:  1 apply filter, 2 use default values, 0 not
%                                  filter.Value:  vector with parameters [low, high, smooth ...]
%                                  filter.Method: i.e. 'circ', 'quadr', 'bandpass'
%                                  filter.Space: 'real' or 'fourier'
%                                  filter.Times: 'apply filter
%                                  n-timesimage = tom_bandpass(im,low,hi,smooth)
%  
%  OUTPUT
%   im_out              filtered image
%
%EXAMPLE
%    im_in=tom_emread('/fs/bmsan/apps/tom_dev/data/2d/mona.em');im_in=im_in.Value;      
%    figure; subplot(1,2,1);tom_imagesc(im_in);
%
%    filter_st.Apply=1;
%    filter_st.Value=[2 4 3];
%    filter_st.Method='circ'; 
%    filter_st.Space='fourier';
%    filter_st.Times=3;                                           
%    imf=tom_apply_filter(im_in,filter_st);   
%    subplot(1,2,2);tom_imagesc(imf);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 25/01/06
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





%select the method

if (filter_st.Apply==0)
    im_out=im_in;
    return;
end;

if (isstruct(im_in));
    error('Only Values accepted ... no Header');
end;


sz=size(im_in);

switch lower(filter_st.Type)

    case 'kernel'
        if (filter_st.Apply==2)
            %Default Settings
            filter_st.Value.radius=round(min(sz)./100);
            filter_st.Value.method='circ';
            filter_st.Value.space='fourier';
            filter_st.Value.times=1;
        end;
            
        for i=1:filter_st.Value.times
            im_out = tom_filter(im_in,filter_st.Value.radius,filter_st.Value.method,filter_st.Value.space);
            im_in=im_out;
        end;
        return;
    
    case 'bandpass'
        if (filter_st.Apply==2)
            %Default Settings
            filter_st.Value.low=0;
            filter_st.Value.high=round(min(sz)./2.2);
            filter_st.Value.smooth=round(min(sz)./15);
            filter_st.Value.times=1;
        end;
        
        %check for empty Values
        if (isempty(filter_st.Value.smooth))
            filter_st.Value.smooth=0;
        end;
        
        for i=1:filter_st.Value.times
            im_out = tom_bandpass(im_in,filter_st.Value.low,filter_st.Value.high,filter_st.Value.smooth);
            im_in=im_out;
        end;
        return;
    
    case 'anisotropic_diffusion' 
        
        
        
    otherwise
            error([filter_st.Type ' not known:']); % bug fixed, SN
            

end;



















% 
% 
% %%% HACK
% try;
% if (strcmp(filter_st.Type,'bandpass')==1)
%     filter_st.Method='bandpass';
%  end;
% end
% %%%%
% 
% 
% if (filter_st.Apply==0)
%     im_out=im_in;
%     return;
% end;
% 
% 
% if (isfield(im_in,'Value'))
%     im_in=im_in.Value;
% end;
% 
% if (filter_st.Apply==2)
%     %calculate default values
%     filter_st.Method='bandpass';
%     filter_st.Value(1)=0;
%     filter_st.Value(2)=round((size(im_in,2)./2));
%     filter_st.Value(3)=0;
%     filter_st.Times=1;
% end;
% 
% if (isfield(filter_st,'Times')==0)
%    filter_st.Times=1;
% end;
% 
% 
% if (isempty(filter_st.Times)==1)
%    filter_st.Times=1; 
% end;
% 
% for i=1:filter_st.Times
% 
%     if (strcmp(filter_st.Method,'bandpass')==1)
%         im_out=tom_bandpass(im_in,filter_st.Value(1),filter_st.Value(2),filter_st.Value(3));
%     end;
% 
%     if (strcmp(filter_st.Method,'quadr')==1 | strcmp(filter_st.Method,'circ')==1 )
%         im_out=tom_filter(im_in,filter_st.Value(1),filter_st.Method,filter_st.Space);
%     end;
% 
%     im_in=im_out;
%     
% end;
