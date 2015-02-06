function output = tom_mark_applyTransf(mode, input, transf, varargin)
%TOM_MARK_APPLYTRANSF applies transformations on a markerset or an image.
% 
%   output = TOM_MARK_APPLYTRANSF(mode, input, transf,
%                                        imreadbinning, varargin)
%   
%   This funciton applies different types of transformations on a image.
%   This can be a specific filter or funcitons, which change the
%   imagedimension and orientiation too.
%   Appart from applying the transformation on the image, it can also be
%   applied on a markerset or on the params for matching. For example if
%   you scale the input image, also the coordinates of the markerset have
%   to be multiplied with a scalar, so that they still match to the
%   transformed image.
%   How the "input" is interpreted, depends on the parameter "mode"
%
%   BE CONSCIOUS, that the em-filetype, as used by the tom toolbox, holds
%   an image transposed. You see this for example, when you display the
%   image using "imagesc". Then you have to handle im.Value' to the function,
%   instead, of im.Value (tom_imagesc  transposes the image implizit).
%   As long as the transformation does not depend on the alignment of the
%   inputimage (for example common filters or resizing with preserving the aspect
%   ratio) it does not matter how the image is passed to
%   TOM_MARK_APPLYTRANSF. Otherwise be carefull!
%
%   The modes which interpret the "input" as markerfile can take:
%     - A filename of an em-markerfile.
%     - A structure with the markerpositons saved in the field "Value".
%     - The array with coordinates as given in 'Value'.
%   If the markerset contains two rows the first row is interpreted as x
%   coordinate and the second is interpreted as y coordinate. If it has
%   more rows, the second and third row are taken as [x; y].
%   The markerset which is returned has the same format as the input. Beware
%   that only those positions in the markerset change, where both the x and
%   y coordinate are greater than 0.
%
%PARAMETERS
%  INPUT
%     mode: selects what type the input is and what will be trasformed.
%         Currently the following modes are supported:
%           - image_pre: Input must be an image (2d-array).
%           - markerset_pre: markerset will be transformed so that the
%             returned markercoordinates match the image transfromed with
%             image_pre.
%           - matchparam_pre: A structure as used in TOM_MARK_FINDMATCHES.
%           - markerset_post: The inverse from markerset_pre to transform
%             the markerset to match the original image.
%     input: The source (e.g. the image or the markerset).
%     transf: Is an array of the structure. Each element has the fields
%         'Apply', 'Value' and 'Type'; in accordance to
%         TOM_REC3DSETFILTERGUI.
%         The transformations are applied in the same order as they appear
%         in the array.
%     further parameters: Depending on the mode the following additional
%         parameters are recognized: 
%         - image_pre: NONE
%         - markerset_pre: imreadbinning = varargin{1}.
%         - matchparam_pre: imreadbinning = varargin{1}.
%         - markerset_post: imreadbinning = varargin{1}.
%         imreadbinning: Specifies that the image was binned already while
%           reading with tom_emreadc.
%  OUTPUT:
%     output: The transformed "input".
%
%   See also TOM_REC3DSETFILTERGUI, TOM_MARK_FINDMATCHES
%
%REFERENCES
%
%   20/03/07
%
%   created by Thomas Haller
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
% %   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch (mode) 
    case {'image_pre'}
        % The input is an image, which has to be transformed.
        % Additional parameters are: NONE
        if (~isempty(varargin))
            error('mode "image_pre" does not allow additional parameters.');
        end;
        output = input;
        for (i=1:length(transf))
            switch (transf(i).Type)
                case {'bandpass', 'kernel'}
                    output = tom_apply_filter(output, transf(i));
                case {'imresize'}
                    if (transf(i).Apply == 1)
                        output = imresize(output, transf(i).Value.scale);
                    end;
                case {'fhandle'}
                    if (transf(i).Apply == 1)
                        vars = transf(i).Value.image_pre_varargin;
                        vars{transf(i).Value.image_pre_vararginpos} = output;
                        [output_tmp{1:nargout(transf(i).Value.image_pre_hfct)}] = feval(transf(i).Value.image_pre_hfct, vars{:});
                        output = output_tmp{transf(i).Value.image_pre_varargoutpos};
                    end;                    
                case {'tom_bin'}
                    if (transf(i).Apply == 1)
                        output = tom_bin(output, transf(i).Value.nbin);
                    end;
                otherwise
                    error(['unrecognized transformation of type "' fransf(i).Type '"']);
            end;            
        end;

    case {'markerset_pre'}
        % The input is a original markerset, which has to be transformed to
        % match an image which was read with binning imreadbinning and which 
        % is transformed with "image_pre"
        % Additional parameters are: imreadbinning.
        if (~isempty(varargin))
            if (length(varargin) > 1)
                error(['The mode markerset_pre only allows one extra parameter (imreadbinning).']);
            end;
            imreadbinning = varargin{1};
            if (~isnumeric(imreadbinning) || numel(imreadbinning)~=1 || uint8(imreadbinning)~=imreadbinning)
                error('The parameter "imreadbinning" must be an integer.');
            end;
        else 
            imreadbinning = 0;
        end;
        
        if (ischar(input))
            input = tom_emread(input);
            is_emstr = true;
        else
            is_emstr = isstruct(input);
        end;
        
        if (is_emstr)
            Value = input.Value;
        else
            Value = input;
        end;
        msize1 = size(Value, 1);
        if (msize1 ~= 2)
            markerset = Value([2 3], :, :);
        else
            markerset = Value;
        end;
        
        idx_existing_marker = all(markerset > 0, 1);
        
        
        markerset(:,idx_existing_marker) = markerset(:,idx_existing_marker) / pow2(imreadbinning);
        for (i=1:length(transf))
            switch (transf(i).Type)
                case {'bandpass', 'kernel'}
                    % do nothing
                case {'imresize'}
                    if (transf(i).Apply == 1)
                        markerset(:,idx_existing_marker) = markerset(:,idx_existing_marker) * transf(i).Value.scale;
                    end;
                case {'fhandle'}
                    if (transf(i).Apply == 1)
                        vars = transf(i).Value.markerset_pre_varargin;
                        vars{transf(i).Value.markerset_pre_vararginpos} = markerset(:,idx_existing_marker);
                        [output_tmp{1:nargout(transf(i).Value.markerset_pre_hfct)}] = feval(transf(i).Value.markerset_pre_hfct, vars{:});
                        markerset(:,idx_existing_marker) = output_tmp{transf(i).Value.markerset_pre_varargoutpos};
                    end;                    
                case {'tom_bin'}
                    if (transf(i).Apply == 1)
                        markerset(:,idx_existing_marker) = markerset(:,idx_existing_marker) / pow2(transf(i).Value.nbin);
                    end;
                otherwise
                    error(['unrecognized transformation of type "' transf(i).Type '"']);
            end;            
        end;
        
        if (msize1 ~= 2)
            output = Value;
            output([2, 3],:,:) = markerset;
        else
            output = markerset;
        end;

        if (is_emstr)
            Value2 = output;
            output = input;
            output.Value = Value2;
        end;
    
    case {'markerset_post'}
        
        if (~isempty(varargin))
            if (length(varargin) > 1)
                error(['The mode markerset_post only allows one extra parameter (imreadbinning).']);
            end;
            imreadbinning = varargin{1};
            if (~isnumeric(imreadbinning) || numel(imreadbinning)~=1 || uint8(imreadbinning)~=imreadbinning)
                error('The parameter "imreadbinning" must be an integer.');
            end;
        else 
            imreadbinning = 0;
        end;
        
        if (ischar(input))
            input = tom_emread(input);
            is_emstr = true;
        else
            is_emstr = isstruct(input);
        end;
        
        if (is_emstr)
            Value = input.Value;
        else
            Value = input;
        end;
        msize1 = size(Value, 1);
        if (msize1 ~= 2)
            markerset = Value([2 3], :, :);
        else
            markerset = Value;
        end;
        
        idx_existing_marker = all(markerset > 0, 1);
        
        
        for (i=1:length(transf))
            switch (transf(i).Type)
                case {'bandpass', 'kernel'}
                    % do nothing
                case {'imresize'}
                    if (transf(i).Apply == 1)
                        markerset(:,idx_existing_marker) = markerset(:,idx_existing_marker) / transf(i).Value.scale;
                    end;
                case {'fhandle'}
                    if (transf(i).Apply == 1)
                        vars = transf(i).Value.markerset_post_varargin;
                        vars{transf(i).Value.markerset_post_vararginpos} = markerset(:,idx_existing_marker);
                        [output_tmp{1:nargout(transf(i).Value.markerset_post_hfct)}] = feval(transf(i).Value.markerset_post_hfct, vars{:});
                        markerset(:,idx_existing_marker) = output_tmp{transf(i).Value.markerset_post_varargoutpos};
                    end;                    
                case {'tom_bin'}
                    if (transf(i).Apply == 1)
                        markerset(:,idx_existing_marker) = markerset(:,idx_existing_marker) * pow2(transf(i).Value.nbin);
                    end;
                otherwise
                    error(['unrecognized transformation of type "' transf(i).Type '"']);
            end;            
        end;
        markerset(:,idx_existing_marker) = markerset(:,idx_existing_marker) * pow2(imreadbinning);
        
        if (msize1 ~= 2)
            output = Value;
            output([2, 3],:,:) = markerset;
        else
            output = markerset;
        end;

        if (is_emstr)
            Value2 = output;
            output = input;
            output.Value = Value2;
        end;
    case {'matchparam_pre'}
        if (~isempty(varargin))
            if (length(varargin) > 1)
                error(['The mode matchparam_pre only allows one extra parameter (imreadbinning).']);
            end;
            imreadbinning = varargin{1};
            if (~isnumeric(imreadbinning) || numel(imreadbinning)~=1 || uint8(imreadbinning)~=imreadbinning)
                error('The parameter "imreadbinning" must be an integer.');
            end;
        else 
            imreadbinning = 0;
        end;
        
        input_w_is_pixeldependent = isnumeric(input.w) && numel(input.w)==1 && input.w > 1;
        input_max_disparity_is_pixeldependent = isnumeric(input.max_disparity) && numel(input.max_disparity)==1 && input.max_disparity > 0;
        input_im_shift_is_pixeldependent = isnumeric(input.im_shift) && ~isempty(input.im_shift);
        
        if (input_w_is_pixeldependent)
            input.w = input.w / pow2(imreadbinning);
        end;
        if (input_max_disparity_is_pixeldependent)
            input.max_disparity = input.max_disparity / pow2(imreadbinning);
        end;
        if (input_im_shift_is_pixeldependent)
            input.im_shift = input.im_shift / pow2(imreadbinning);
        end;  
        
        for (i=1:length(transf))
            switch (transf(i).Type)
                case {'bandpass', 'kernel'}
                    % do nothing
                case {'imresize'}
                    if (transf(i).Apply == 1)
                        if (input_w_is_pixeldependent)
                            input.w = input.w * transf(i).Value.scale;
                        end;
                        if (input_max_disparity_is_pixeldependent)
                            input.max_disparity = input.max_disparity * transf(i).Value.scale;
                        end;
                        if (input_im_shift_is_pixeldependent)
                            input.im_shift = input.im_shift * transf(i).Value.scale;
                        end;  
                    end;
                case {'fhandle'}
                    if (transf(i).Apply == 1)
                        inputfieldnames = fieldnames(input);
                        existing_field = true(1,length(transf(i).Value.matchparam_pre_fields));
                        for (ifield=1:length(inputfieldnames))
                            fieldname = inputfieldnames{ifield};
                            jfield = find(strcmp(fieldname, transf(i).Value.matchparam_pre_fields));
                            if (~isempty(jfield))
                                existing_field(jfield) = false;
                                if (length(jfield) > 1)
                                    warning(['The transformation-structure at index #' num2str(i) ' contains more times the value "' fieldname '"']);
                                    jfield = jfield(1);
                                end;
                                
                                if (~exist(['input_' fieldname '_is_pixeldependent'], 'var'))
                                    warning(['The transformation-structure at index #' num2str(i) ' contains a reference to the field "' fieldname '" but this parameter does not depend on the actual pixelsize']);
                                elseif (eval(['input_' fieldname '_is_pixeldependent']))
                                    vars = transf(i).Value.matchparam_pre_varargin(jfield, :);
                                    vars{transf(i).Value.matchparam_pre_vararginpos(jfield)} = input.(fieldname);
                                    hfct = transf(i).Value.matchparam_pre_hfct{jfield};
                                    [output_tmp{1:nargout(hfct)}] = feval(hfct, vars{:});
                                    input.(fieldname) = output_tmp{transf(i).Value.matchparam_pre_varargoutpos(jfield)};
                                end;
                            else
                                if (exist(['input_' fieldname '_is_pixeldependent'], 'var'))
                                    warning(['The matchparam-field ' fieldname ' depends on the image-transformations, but it has no quivalent in the transformation-structure.']);
                                end;
                            end;
                        end;
                        if (any(existing_field))
                            warning(['The transformation-structure at index #' num2str(i) ' contains fieldnames which do not occure in the matchparam structure.']);
                        end;
                        
                    end;                    
                case {'tom_bin'}
                    if (transf(i).Apply == 1)
                        if (input_w_is_pixeldependent)
                            input.w = input.w / pow2(transf(i).Value.nbin);
                        end;
                        if (input_max_disparity_is_pixeldependent)
                            input.max_disparity = input.max_disparity / pow2(transf(i).Value.nbin);
                        end;
                        if (input_im_shift_is_pixeldependent)
                            input.im_shift = input.im_shift / pow2(transf(i).Value.nbin);
                        end;  
                    end;
                otherwise
                    error(['unrecognized transformation of type "' transf(i).Type '"']);
            end;            
        end;        
        if (input_w_is_pixeldependent && input.w <= 1)
            input.w = 1.1;
        end;
        output = input;
        
    
    otherwise
        error(['option "' mode '" is not recognized']);    
end;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function transf = example_call()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transf(1).Type = 'fhandle';
transf(end).Apply = 1;
transf(end).Value.image_pre_varargin = {[], 0.5, 'cubic'};
transf(end).Value.image_pre_vararginpos = 1;
transf(end).Value.image_pre_varargoutpos = 1;
transf(end).Value.image_pre_hfct = @imresize;
transf(end).Value.markerset_pre_varargin = {[], transf(1).Value.image_pre_varargin{2}};
transf(end).Value.markerset_pre_vararginpos = 1;
transf(end).Value.markerset_pre_varargoutpos = 1;
transf(end).Value.markerset_pre_hfct = @times;
transf(end).Value.matchparam_pre_fields = {'w', 'im_shift', 'max_disparity'};
transf(end).Value.matchparam_pre_varargin = {[], transf(1).Value.image_pre_varargin{2}; ...
                                             [], transf(1).Value.image_pre_varargin{2}; ...
                                             [], transf(1).Value.image_pre_varargin{2}; ...
                                             [], transf(1).Value.image_pre_varargin{2}; ...
                                             [], transf(1).Value.image_pre_varargin{2}; };
transf(end).Value.matchparam_pre_vararginpos = [1 1 1 1 1];
transf(end).Value.matchparam_pre_varargoutpos = [1 1 1 1 1];
transf(end).Value.matchparam_pre_hfct = {@times, @times, @times,@times, @times};
transf(end).Value.markerset_post_varargin = {[], 1/transf(1).Value.image_pre_varargin{2}};
transf(end).Value.markerset_post_vararginpos = 1;
transf(end).Value.markerset_post_varargoutpos = 1;
transf(end).Value.markerset_post_hfct = @times;


transf(2).Type = 'tom_bin';
transf(end).Apply = 1;
transf(end).Value.nbin = 1;

transf(3).Type = 'bandpass';
transf(end).Apply = 2;
transf(end).Value.times = 1;
transf(end).Value.low = 10;
transf(end).Value.high = 100;
transf(end).Value.smooth = 0;
transf(end).Value.space = 'real';
transf(end).Value.method = 'quadr';
transf(end).Value.radius = 0;
 
transf(end+1).Type = 'imresize';
transf(end).Apply = 1;
transf(end).Value.scale = 2;
