classdef overlayImg < handle
    % overlays two datasets and displays them as overlayed images
    % call examples: 
    % c = overlayImg(backimg,overimg)
    % c = overlayImg(backimg,overimg,[fromRange toRange],...
    %                        [fromTransparent toTransparent])
    %
    % after the class is constructed you can edit the properties either
    % directly and call update() manually or you can use the setPROPERTY
    % functions and the figure will update automatically. 
    %
    % It is possible to set the value ranges, the color maps, a
    % transparency range and the alpha value.
    %
    % Where backimg is the image in the back, overimg the image overlay.
    % If the third parameter is given the values of overimg are limited to
    % the range given. If the fourth parameter is given the values in the
    % given range are transparent (for example for just showing negative and
    % positive values, but not around zero).
    %
    % If you want to display indexed images set the corresponding range 
    % (backRange or overRange) to [ 1 length(cmap) ], cmap being the
    % colormap for the image.
    % 
    % The datatip will show coordinates, the value of the background and
    % overlay image.
    %
    % usage example:
    %   rgb = imread('ngc6543a.jpg');  % reading an image
    %   grayimg = mean(double(rgb),3); % make a grayscale image
    %   s = size(grayimg);
    %   [x,y]=meshgrid(1:s(2),1:s(1));  % make grid    
    %   oimg = sin(y/100).^2+cos(x/250).^2;   % some meaningless overlay
    %   c = overlayImg(grayimg,oimg);  % overlay images
    %
    %   % from now on we can interact with the object
    %   % let's check the values in the background and overlay image
    %   c.showColorBars()
    %   % now setting up a transparent range
    %   c.setTranspRange([0 0.5])
    %   % reduce range in background image
    %   c.setBackRange([0 150])
    %   % adjust alpha value
    %   c.setAlpha(0.5)
    %   % don't show overlay values above 1.5
    %   c.setOverRange([0 1.5])
    
    % Jochen Deibele 2009, 2010
    
    properties
        mFigHdl = [];           % figure handle
        mCbHdlb = [];           % background colorbar handle
        mCbHdlo = [];           % overlay colorbar handle
        
        % background image
        backImg = [];           % background image
        backMap = gray(256);    % color map background image
        backRange = [];         % background image range = min:max
        
        % overlay
        overImg = [];           % overlay image
        %overMap = jet(256);     % Color map overlay
        overMap = gray;
        overRange = [];         % imaging range empty = min:max
        
        % transparancy 
        transpRange = [];       % [from to] from < overImg < to will not be shown
        alpha = 0.7;            % alpha value
    end
    methods (Static)
        function rRange = detRange(aIm, aRange)
            % aRange = detRange(aIm, aRange)
            % convenience function to determine value range. If aRange is
            % empty the range will be taken from min(aIm(:)) and
            % max(aIm(:)) ignoring values equal to +-Inf.
            if isempty(aRange)
                aIm_mod = aIm;
                aIm_mod(isinf(aIm)) = [];
                rRange = [min(aIm_mod(:)) max(aIm_mod(:))];
            else
                rRange = aRange;
            end;
        end;
            
        function si = idxImg(aIm, aRange, aMap)
            % si = idxImg(aIm, aRange, aMap)
            % converts a (not indexed) image to an indexed image in the
            % given colormap, saturating if values exceed aRange
            aIm = double(aIm);
            aRange = overlayImg.detRange(aIm, aRange);  % determine range
            si = (aIm - aRange(1)) / (aRange(2) - aRange(1)) * (size(aMap,1)-1)+1;
            si = floor(si);
            si(aIm < aRange(1)) = 1;
            si(aIm > aRange(2)) = size(aMap,1);
        end;
    end;
    methods
        function obj = overlayImg(aBackImg, aOverImg, aOverRange, aTranspRange) 
            % overlayImg(aBackImg, aOverImg, aOverRange, aTranspRange) 
            % constructor
            % aBackImg : background image
            % aOverImg : overlay image
            % aOverRange (optional): Limit the range of the overlay image
            % aTranspRange: Range of values for the overlay image which is
            %               not shown
            
            % initialize values 
            obj.backImg = aBackImg;
            obj.overImg = aOverImg;
            if nargin > 2
                obj.overRange = aOverRange;
            end;
            if nargin > 3
                obj.transpRange = aTranspRange;
            end;
            obj.mFigHdl = figure();
            axes();
            
            % change data tip callback to own function
            dcm_obj = datacursormode(obj.mFigHdl);
            set(dcm_obj,'UpdateFcn',@obj.datatipCB);
            
            obj.update();
        end
        
        function tiptext = datatipCB(obj, cbobj, event_obj)
            % datatip callback function
            % the datatip shows coordinates and the value of the background
            % as well as of the overlay image
            posx = round(event_obj.Position(1));
            posy = round(event_obj.Position(2));
            tiptext = sprintf('x: %i y: %i\nback: %2.2f\nover: %2.2f',posx,posy,obj.backImg(posy,posx),obj.overImg(posy,posx));
        end;
        
        function update(obj)
            % renders the image
            
            % create common colormap
            map = [obj.backMap;obj.overMap];
            
            % convert to indexed image
            gi = obj.idxImg(obj.backImg, obj.backRange, obj.backMap);
            oi = obj.idxImg(obj.overImg, obj.overRange, obj.overMap);
            
            ax = get(obj.mFigHdl,'Children');
            image(gi,'Parent',ax)
            hold(ax,'all');
            h = image(oi + size(obj.backMap,1),'Parent',ax);
            hold(ax,'off');
            colormap(ax,map);
            
            % is a transparent range given?
            if isempty(obj.transpRange)
                % no, we use the same alpha value for all of the overlay
                % image
                set(h,'AlphaData',obj.alpha);
            else
                % yes, we have to make the range transparent
                alMask = zeros(size(obj.overImg));
                alMask(obj.overImg<double(obj.transpRange(1))) = 1;
                alMask(obj.overImg>double(obj.transpRange(2))) = 1;
                set(h,'AlphaData',alMask*obj.alpha);
            end;
            
            % if colorbars are shown update them
            if ~isempty(obj.mCbHdlb)
                obj.showColorBars();
            end;
        end
       
        function showColorBars(obj)
            % background colorbar
            aHdl = obj.mCbHdlb;
            % check if colorbar is already on screen, if not create new
            if isempty(aHdl) || ~ishandle(aHdl)
                aHdl = figure;
                pos = get(aHdl,'Position');
                set(aHdl,'Position',[pos(1)-pos(3)*2/5-15 pos(2) pos(3)/5 pos(4)]);
            else
                figure(aHdl);
            end;
            obj.mCbHdlb = aHdl;
            % determine value range 
            aRange = overlayImg.detRange(obj.backImg, obj.backRange);
            aMap = obj.backMap;
            aRangeVec = (aRange(2)-aRange(1))/(size(aMap,1)-1)*(0:(size(aMap,1)-1))+aRange(1);
            imagesc(1,aRangeVec,aRangeVec.');
            set(gca,'YDir','normal');
            colormap(aMap);
            title('legend back');
            
            % overlay colorbar 
            aHdl = obj.mCbHdlo;
            % check if colorbar is already on screen, if not create new
            if isempty(aHdl) || ~ishandle(aHdl)
                aHdl = figure;
                pos = get(aHdl,'Position');
                set(aHdl,'Position',[pos(1)-pos(3)/5-7 pos(2) pos(3)/5 pos(4)]);
            else
                figure(aHdl);
            end;
            obj.mCbHdlo = aHdl;
            
            % determine value range 
            aRange = overlayImg.detRange(obj.overImg, obj.overRange);            
            aMap = obj.overMap;
            aRangeVec = (aRange(2)-aRange(1))/(size(aMap,1)-1)*(0:(size(aMap,1)-1))+aRange(1);
            imagesc(1,aRangeVec,aRangeVec.');
            set(gca,'YDir','normal');
            colormap(aMap);
            title('legend over');
            % focus back to figure
            figure(obj.mFigHdl);
        end;
        function setBackImg(obj, aValue)
            obj.backImg = aValue;
            obj.update();
        end
        function setOverImg(obj, aValue)
            obj.overImg = aValue;
            obj.update();
        end
        function setAlpha(obj, aValue)
            obj.alpha = aValue;
            obj.update();
        end;
        function setOverRange(obj, aValue)
            obj.overRange = aValue;
            obj.update();
        end;
        function setBackRange(obj, aValue)
            obj.backRange = aValue;
            obj.update();
        end;
        function setOverMap(obj, aValue)
            obj.overMap = aValue;
            obj.update();
        end;
        function setBackMap(obj, aValue)
            obj.backMap = aValue;
            obj.update();
        end;
        function setTranspRange(obj, aValue)
            obj.transpRange = aValue;
            obj.update();
        end;
    end
end