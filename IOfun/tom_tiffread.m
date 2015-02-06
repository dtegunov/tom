function Data = tom_tiffread( tiff_name )
%TOM_TIFFREAD reads data in TIFF-file format
%
%   Data = tom_tiffread(filename);
%
%   Reads a tiff-Image File and the TVIPS specific Header.
%   The Header Information makes it possible to read
%   unsigned and signed tiff-Image Files.
%   The TVIPS Header Information is coded in the registerd
%   private Tiff-Tags 37706 - 37710.
%   Details on the Header Information can be viewed in
%   the TVIPS TemData Extension documentation(tiff.html).
%
%  INPUT
%   filename
%
%  OUTPUT
%   Data            Structure of Image Data
%   Data.Value      Image Data
%   Data.Header     Header Information
%
%EXAMPLE
%               i=tom_tiffread;
%                   a fileselect-box appears and the TIFF-file can be picked
%
%               i=tom_tiffread('Image.tif');
%


error(nargchk(0,1,nargin))
if nargin <1 
    [filename, pathname] = uigetfile({'*.tif';'*.*'}, 'Pick a TIFF-file');
    if isequal(filename,0) | isequal(pathname,0) disp('No data loaded.'); return; end;
    tiff_name = [pathname filename];
end;

[Image, Header] = tom_tiffreadmx(tiff_name);
Data=struct('Value', Image, 'Header', Header);