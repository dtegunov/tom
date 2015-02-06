function tom_av2_check_htl_file(text_file,increment,filter)
%TOM_AV2_CHECK_HTL_FILE(filename)  manual ispection of htl files
%
%   tom_av2_check_htl_file(filename)
%
%   TOM_AV2_CHECK_HTL_FILE manual ispection of htl files
%   
%  
%
%PARAMETERS
%
%  INPUT
%   text_file           htl fiel
%   increment           (500) check every n-th file
%   filter              apply a kernel filter
%
%  OUTPUT
%  
%
%
%EXAMPLE
%   
% tom_av2_check_htl_file('16_corrf_high_128.mat_16_corrf_low_128.mat.htl',1 00,2);
%  
%
% 
%
%
%

if (nargin <3)
    filter=1;
end;

if (nargin <2)
    increment=400;
end;


d=importdata(text_file);
figure;

try
    ll=length(d.textdata);
catch
    ll=length(d);
end;


for i=1:increment:ll
    try
        p1=d.textdata{i,1};
        p2=d.textdata{i,2};
    catch ME
        [a b]=strtok(d{i},' ');
        p1=strrep(a,' ','');
        p2=strrep(strrep(b,' ',''),',','');
    end;
    
    
    im1=tom_spiderread(p1);
    im2=tom_spiderread(p2);
    subplot(1,2,1); tom_imagesc(tom_filter(im1.Value,filter));
    subplot(1,2,2); tom_imagesc(tom_filter(im2.Value,filter));
    ginput(1);
end;

