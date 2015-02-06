% Used to convert an em-structure (with the fields 'Value' and 'Header'
% into the bits and bits of the em-header. Since the is not always
% consistent, the function tries to figure out the correct meaning of each
% header field.
function [magic, comment, emdata, userdata] = tom_emheader_struct_2_data(st, data_type)


if (isnumeric(st))
    if (isempty(st))
        st = struct('Header', {struct()});
    else 
        st = struct('Value', {st}, 'Header', {struct()});
    end;
elseif (isstruct(st))
    if (numel(st)==1 && isfield(st,'Value') && isnumeric(st.Value) && isfield(st,'Header') && isstruct(st.Header) && numel(st.Header)==1)
    else
        st = struct('Header', {st});
    end;
else
    error('st is not valid EM-data');
end;

if (exist('data_type','var') && ~isempty(data_type) && ~(ischar(data_type)&&isempty(data_type)))
    if (~ischar(data_type))
        error('Datatype must be a string (int8, int16, int32, single, complex32, double, complex64)');
    end;
    switch lower(data_type)
        case 'int8'
            data_type = 1;
        case 'int16'
            data_type = 2;
        case 'int32'
            data_type = 4;
        case {'single', 'float', 'short'}
            data_type = 5;
        case {'complex32', 'single_complex', 'float_complex'}
            data_type = 8;
        case 'double'
            data_type = 9;
        case {'complex64', 'double_complex'}
            data_type = 10;
        otherwise
            error('Unexpected datatype (only int8, int16, int32, single, complex32, double, complex64 allowed)');
    end;
else
    data_type = [];
    if (isfield(st, 'Value'))
        if (any(imag(st.Value)))
            switch (class(st.Value))
                case 'single'
                    data_type = 8;
                case 'double'
                    data_type = 10;
                otherwise
                    error('Complex type only allowed for single.');
            end;
        else
            switch (class(st.Value))
                case 'int8'
                    data_type = 1;
                case 'int16'
                    data_type = 2;
                case 'int32'
                    data_type = 4;
                case 'single'
                    data_type = 5;
                case 'double'
                    data_type = 9;
                otherwise
                    error('Data-type not available in EM.');
            end;
        end;        
    end;
end;




Header = st.Header;



magic = int8([6, 0, 0, 5]); % DEFAULT: PC, SINGLE
if (isfield(Header, 'EM') && isfield(Header.EM, 'Magic') && numel(Header.EM.Magic)==4)
    magic(1:4) = int8(Header.EM.Magic(1:4));
end;
if (~isempty(data_type))
    magic(4) = data_type;
end;



if (isfield(Header, 'EM'))
    if (isfield(Header.EM, 'Comment'  ) && isnumeric(Header.EM.Comment  ) && numel(Header.EM.Comment  )== 80) comment  = int8 (Header.EM.Comment  (1:80));  end;
    if (isfield(Header.EM, 'Parameter') && isnumeric(Header.EM.Parameter) && numel(Header.EM.Parameter)== 40) emdata   = int32(Header.EM.Parameter(1:40));  end;
    if (isfield(Header.EM, 'Fillup'   ) && isnumeric(Header.EM.Fillup   ) && numel(Header.EM.Fillup   )==256) userdata = int8 (Header.EM.Fillup   (1:256)); end;
end;
if (~exist('comment' , 'var')) comment  = zeros(1, 80, 'int8' );          end;
if (~exist('emdata'  , 'var')) emdata   = zeros(1, 40, 'int32');          end;
if (~exist('userdata', 'var')) userdata = zeros(1,256, 'int8' );          end;

if (isfield(Header, 'Voltage'))             emdata(1) = Header.Voltage;                             end;
if (isfield(Header, 'Cs'))                  emdata(2) = Header.Cs * 1000;                           end;
if (isfield(Header, 'Aperture'))            emdata(3) = Header.Aperture;                            end;
if (isfield(Header, 'Magnification'))       emdata(4) = Header.Magnification;                       end;
if (isfield(Header, 'Postmagnification'))   emdata(5) = Header.Postmagnification * 1000;            end;
if (isfield(Header, 'Exposuretime'))        emdata(6) = Header.Exposuretime * 1000;                 end;
if (isfield(Header, 'Objectpixelsize'))     emdata(7) = Header.Objectpixelsize * 1000;              end;
if (isfield(Header, 'Microscope'))          emdata(8) = tom_em_microscope_name(Header.Microscope);  end;
if (isfield(Header, 'Pixelsize'))           emdata(9) = Header.Pixelsize * 1000;                    end;
if (isfield(Header, 'CCDArea'))             emdata(10) = Header.CCDArea * 1000;                     end;
if (isfield(Header, 'Defocus'))             emdata(11) = Header.Defocus;                            end;
if (isfield(Header, 'Astigmatism'))         emdata(12) = Header.Astigmatism;                        end;
if (isfield(Header, 'AstigmatismAngle'))    emdata(13) = Header.AstigmatismAngle * 1000;            end;
if (isfield(Header, 'FocusIncrement'))      emdata(14) = Header.FocusIncrement;                     end;
if (isfield(Header, 'CountsPerElectron'))   emdata(15) = Header.CountsPerElectron * 1000;           end;
if (isfield(Header, 'Intensity'))           emdata(16) = Header.Intensity * 1000;                   end;
if (isfield(Header, 'EnergySlitwidth'))     emdata(17) = Header.EnergySlitwidth;                    end;
if (isfield(Header, 'EnergyOffset'))        emdata(18) = Header.EnergyOffset;                       end;
if (isfield(Header, 'Tiltangle'))           emdata(19) = Header.Tiltangle * 1000;                   end;
if (isfield(Header, 'Tiltaxis'))            emdata(20) = Header.Tiltaxis * 1000;                    end;
if (isfield(Header, 'Marker_X'))            emdata(23) = Header.Marker_X;                           end;
if (isfield(Header, 'Marker_Y'))            emdata(24) = Header.Marker_Y;                           end;
if (isfield(Header, 'Username'))
    u = Header.Username(1:min(length(Header.Username),20));
    userdata( 1:20) = [u(:); zeros(20-length(u),1)];
end;
if (isfield(Header, 'Date'))
    u = Header.Date(1:min(length(Header.Date),8));
    userdata(21:28) = [u(:); zeros(8-length(u),1)];
end;


