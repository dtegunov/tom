function varargout = tom_mark_hist(handle)
%IMCONTRAST Adjust Contrast tool.
%   IMCONTRAST creates an Adjust Contrast tool in a separate figure that is
%   associated with the grayscale image in the current figure, called the 
%   target image. The Adjust Contrast tool is an interactive contrast and 
%   brightness adjustment tool that you can use to adjust the
%   black-to-white mapping used to display the image. The tool works by
%   modifying the CLim property.
%
%   Note: The Adjust Contrast tool can handle grayscale images of class 
%   double and single with data ranges that extend beyond the default
%   display range, which is [0 1]. For these images, IMCONTRAST sets the
%   histogram limits to fit the image data range, with padding at the upper
%   and lower bounds.
%
%   IMCONTRAST(H) creates an Adjust Contrast tool associated with the image
%   specified by the handle H. H can be an image, axes, uipanel, or figure
%   handle. If H is an axes or figure handle, IMCONTRAST uses the first
%   image returned by FINDOBJ(H,'Type','image').
%
%   HFIGURE = IMCONTRAST(...) returns a handle to the Adjust Contrast tool
%   figure.
%
%   Remarks
%   -------
%   The Adjust Contrast tool presents a scaled histogram of pixel values
%   (overly represented pixel values are truncated for clarity). Dragging
%   on the left red bar in the histogram display changes the minimum value.
%   The minimum value (and any value less than the minimum) displays as
%   black. Dragging on the right red bar in the histogram changes the
%   maximum value. The maximum value (and any value greater than the
%   maximum) displays as white. Values in between the red bars display as
%   intermediate shades of gray.
%
%   Together the minimum and maximum values create a "window". Stretching
%   the window reduces contrast. Shrinking the window increases contrast.
%   Changing the center of the window changes the brightness of the image.
%   It is possible to manually enter the minimum, maximum, width, and
%   center values for the window. Changing one value automatically updates
%   the other values and the image.
%
%   Window/Level Interactivity
%   --------------------------
%   Clicking and dragging the mouse within the target image interactively
%   changes the image's window values. Dragging the mouse horizontally from
%   left to right changes the window width (i.e., contrast). Dragging the
%   mouse vertically up and down changes the window center (i.e.,
%   brightness). Holding down the CTRL key when clicking accelerates
%   changes. Holding down the SHIFT key slows the rate of change. Keys must
%   be pressed before clicking and dragging.
%
%   Example
%   -------
%
%       imshow('pout.tif')
%       imcontrast(gca)
%
%    See also IMADJUST, IMTOOL, STRETCHLIM.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.1.8.13 $  $Date: 2006/11/08 17:49:22 $


% Do sanity checking on handles and take care of the zero-argument case.

eid = sprintf('Images:%s:noImageInFigure', mfilename);
if (nargin == 0)
    handle = get(0, 'CurrentFigure');
    if isempty(handle)
        error(eid, '%s expects a current figure containing an image.', ...
            upper(mfilename));
    end
end

iptcheckhandle(handle, {'figure', 'axes', 'image', 'uipanel'}, mfilename, ...
    'H', 1);

[imageHandle, axHandle, figHandle] = imhandles(handle);

if (isempty(imageHandle))
    error(eid, 'The figure must contain at least one image.')
end

% IMCONTRAST uses the first image if there are multiple.
imageHandle = imageHandle(1);
axHandle = axHandle(1);

imgModel = validateImage(imageHandle);

% Install pointer manager in the figure containing the target image.
iptPointerManager(figHandle);

% Display the original image.
figure(figHandle);

% Open a new figure or bring up an existing one?
hFig = getappdata(axHandle, 'imcontrastFig');
if (~isempty(hFig))
    figure(hFig);
    if (nargout > 0)
        varargout{1} = hFig;
    end
    return
end

% The default display range for double images is [0 1].  This default setting
% does not work for double images that are really outside this range; users
% would not even see the draggable window on the histogram (geck 227671).  In
% these cases, we throw a warning and set the display range to include the image
% data range.
badDisplayRange = isDisplayRangeOutsideDataRange(imageHandle,axHandle);
if badDisplayRange
    cdata = get(imageHandle, 'CData');
    imageRange = [double(min(cdata(:))) double(max(cdata(:)))];
    response = displayWarningDialog(get(axHandle,'Clim'), imageRange);
    if strcmpi('OK',response)
        set(axHandle,'Clim', imageRange);
        hHistFig = createContrastTool;
    else
        hHistFig = [];
    end
else
    hHistFig = createContrastTool;
end


    %=====================================
    function hHistFig = createContrastTool

        % Create the histogram palette.
        hHistFig = createHistogramPalette(imageHandle, imgModel);

        % Install pointer manager in the contrast tool figure.
        iptPointerManager(hHistFig);

        iptwindowalign(figHandle, 'left', hHistFig, 'left');
        iptwindowalign(figHandle, 'bottom', hHistFig, 'top');

        % Close the histogram palette if the image goes away.
        reactToImageChangesInFig(imageHandle, @deleteFcn);

        %==================
        function deleteFcn
            if ishandle(hHistFig)
                delete(hHistFig);
            end
        end

        set(hHistFig, 'visible', 'on');
    end

if (nargout > 2)
    error('Images:imcontrast:tooManyOutputs', ...
        'Zero or one output arguments expected.');
elseif (nargout == 1)
    varargout{1} = hHistFig;
end

end % imcontrast

%=============================================================================
function hFig = createHistogramPalette(imageHandle, imgModel)

hImageAx = ancestor(imageHandle, 'axes');
hImageFig = ancestor(imageHandle, 'figure');

isCallerIMTOOL = strncmpi('Image Tool', get(hImageFig,'Name'), 10);

% initializing variables for function scope
cbk_id_cell = {};
isDoubleOrSingleData = false;
undoMenu = [];
redoMenu = [];
undoAllMenu = [];
originalImagePointerBehavior = [];

% variables used for enabling keeping a history of changes
climHistory = [];
currentHistoryIndex = 0;

% Calculate properties of the image to be used throughout this function.
isDoubleOrSingleData = any(strmatch(getClassType(imgModel), {'double','single'}));
newClim = getClim;

% boolean variable used to prevent recursing through the event handler
% and duplicate entries in the history
blockEventHandler = false;

% boolean variable used to indicate if window level operation has started
% so that we know when to save the clim.
startedWindowLevel = false;

defaultFigurePos = [0 0 560 300];
hFig = figure('visible', 'off', ...
    'toolbar', 'none', ...
    'menubar', 'none', ...
    'IntegerHandle', 'off', ...
    'NumberTitle', 'off', ...
    'Name', createFigureName('Adjust Contrast',hImageFig), ...
    'HandleVisibility', 'callback', ...
    'units', 'pixels', ...
    'position', defaultFigurePos);

setappdata(hImageAx, 'imcontrastFig', hFig);

createMenubar;

% create a blank uitoolbar to get docking arrow on the mac as a workaround to
% g222793.
if ismac
    h = uitoolbar(hFig);
end

margin = 5;
hFigFlow = uiflowcontainer('Parent', hFig,...
    'FlowDirection', 'TopDown', ...
    'Margin', margin);

% Create panel that contains window edit boxes, auto scaling, and 
% image range.
[backgroundColor clipPanelAPI windowClipPanelWidth] = ...
    createWindowClipPanel(hFigFlow, imgModel);
editBoxAPI = clipPanelAPI.editBoxAPI;
scalePanelAPI = clipPanelAPI.scalePanelAPI;

% Create HistogramPanel.
hPanelHist = imhistpanel(hFigFlow,imageHandle);
set(hPanelHist, 'Tag','HistAxesPanel');

% Turn off HitTest of the histogram so it doesn't intercept button down
% events - g330176
hHistogram = findobj(hPanelHist, 'type', 'line');
set(hHistogram, 'HitTest', 'off');

% Create Draggable Clim Window on the histogram.
hHistAx = findobj(hPanelHist,'type','axes');
histStruct = getHistogramData(imageHandle);
maxCounts = max(histStruct.counts);
windowAPI = createClimWindowOnAxes(hHistAx,newClim,maxCounts);

% Create Status Panel
[hStatusLabel statusPanelWidth]= createStatusPanel(hFigFlow);

setUpCallbacksOnDraggableWindow;
setUpCallbacksOnWindowWidgets;
setUpCallbacksOnAutoScaling;
set(hFig,'Color', backgroundColor);
setChildColorToMatchParent(hPanelHist, hFig);


% Enable window/leveling through the mouse.
origBtnDwnFcn = get(imageHandle, 'ButtonDownFcn');
[startedWindowLevel winLevelCbkStartId winLevelCbkStopId] = ...
    attachWindowLevelMouseActions;

figureWidth = max([statusPanelWidth windowClipPanelWidth]);
set(hFig,'Position', [0 0 figureWidth 300]);
set(hFig, 'DeleteFcn', @closeHistFig);

reactToPropertyChanges;

updateAllAndSaveInHistory(newClim);

    %===============================
    function reactToPropertyChanges

        % Update the window if the CLIM changes from outside the tool.
        axes_handle = handle(hImageAx);
        clim = axes_handle.findprop('Clim');
        ClimListener = handle.listener(axes_handle, clim, ...
            'PropertyPostSet', @updateTool);
    
        %============================
        function updateTool(src,evt)
            
            if blockEventHandler
                return
            end
            if startedWindowLevel
                updateAll(evt.NewValue);
            else
                updateAllAndSaveInHistory(evt.NewValue);
            end
        end
        
        % Close the histogram palette if the image's cdata or cdatamapping changes.
        image_handle = handle(imageHandle);
        cdata = image_handle.findprop('CData');
        CdataListener = handle.listener(image_handle, cdata,...
            'PropertyPostSet', @(varargin) close(hFig));
        cdatamapping = image_handle.findprop('CDataMapping');
        CdataMappingListener = handle.listener(image_handle, cdatamapping,...
            'PropertyPostSet', @(varargin) close(hFig));
        
        setappdata(hFig, 'TargetListener', ...
            [ClimListener CdataListener CdataMappingListener]);
    end

    %======================================================================
    function [startedWindowLevel,winLevelCbkStartId, winLevelCbkStopId] = ...
            attachWindowLevelMouseActions

        startedWindowLevel = false;
        
        % Attach window/level mouse actions.
        if ~isCallerIMTOOL
            set(imageHandle,...
                'ButtonDownFcn', @(hobj,evt)(windowlevel(imageHandle, hFig)));
        end
    
        winLevelCbkStartId = iptaddcallback(imageHandle,...
            'ButtonDownFcn',@winLevelStarted);
        
        winLevelCbkStopId = iptaddcallback(hImageFig,...
            'WindowButtonUpFcn',@winLevelStopped);
        
        % Change the pointer to window/level when over the image.  Remember
        % the original pointer behavior so we can restore it later in
        % closeHistFig.
        originalImagePointerBehavior = iptGetPointerBehavior(imageHandle);
        enterFcn = @(f,cp) set(f, 'Pointer', 'custom', ...
                                  'PointerShapeCData', getWLPointer);
        iptSetPointerBehavior(imageHandle, enterFcn);
        
        %========================================
        function PointerShapeCData = getWLPointer
            iconRoot = ipticondir;
            cdata = makeToolbarIconFromPNG(fullfile(iconRoot, ...
                                                    'cursor_contrast.png'));
            PointerShapeCData = cdata(:,:,1) + 1;

        end

        %================================
        function winLevelStarted(obj,evt)
            startedWindowLevel = true;
        end

        %================================
        function winLevelStopped(obj,evt)
            startedWindowLevel = false;
        end
    end

    %===================================
    function closeHistFig(obj,eventData)
        
        if blockEventHandler
            return;
        end
        if isappdata(hImageAx, 'imcontrastFig')
            rmappdata(hImageAx, 'imcontrastFig');
        end
        targetListeners = getappdata(hFig, 'TargetListener');
        delete(targetListeners);
        iptremovecallback(imageHandle, ...
            'ButtonDownFcn', winLevelCbkStartId);
        iptremovecallback(hImageFig,...
            'WindowButtonDownFcn', winLevelCbkStopId);
        
        % Restore original image pointer behavior.
        iptSetPointerBehavior(imageHandle, originalImagePointerBehavior);

        if ~isCallerIMTOOL
            set(imageHandle, 'ButtonDownFcn', origBtnDwnFcn);
        end
        
        deleteCursorChangeOverDraggableObjs(cbk_id_cell);
    end

    %=====================
    function createMenubar

        filemenu = uimenu(hFig, ...
            'Label', '&File', ...
            'Tag', 'file menu');
        editmenu = uimenu(hFig, ...
            'Label', '&Edit', ...
            'Tag', 'edit menu');

        if isJavaFigure
            uimenu(hFig, ...
                'Label', '&Window', ...
                'tag', 'winmenu', ...
                'Callback', winmenu('callback'));
        end

        % File menu
        uimenu(filemenu, ...
            'Label', '&Close', ...
            'Accelerator', 'W', ...
            'Callback', @(varargin) close(hFig));

        % Edit menu
        undoMenu = uimenu(editmenu, ...
            'Label', '&Undo', ...
            'Accelerator', 'Z', ...
            'Tag', 'undo menu item', ...
            'Callback', @undoLastChange);
        redoMenu = uimenu(editmenu, ...
            'Label', '&Redo', ...
            'Accelerator', 'Y', ...
            'Tag', 'redo menu item', ...
            'Callback',@redoLastUndo);
        undoAllMenu = uimenu(editmenu, ...
            'Label', 'Undo &All', ...
            'Separator', 'on', ...
            'Tag', 'undo all menu item', ...
            'Callback', @undoAllChanges);

        % Help menu
        if ~isdeployed
            helpmenu = uimenu(hFig, ...
                'Label', '&Help', ...
                'Tag', 'help menu');
            invokeHelp = @(varargin) ...
                ipthelp('imtool_imagecontrast_help','Adjust Contrast'); 
            uimenu(helpmenu, ...
                'Label', 'Adjust Contrast Help', ...
                'Tag', 'imcontrast help menu', ...
                'Callback', invokeHelp);
            iptstandardhelp(helpmenu);
        end
    end % createMenubar

    %=======================================
    function setUpCallbacksOnDraggableWindow

        buttonDownTable = {
            windowAPI.centerLine.handle  @centerPatchDown;
            windowAPI.centerPatch.handle @centerPatchDown;
            windowAPI.maxLine.handle     @minMaxLineDown;
            windowAPI.minLine.handle     @minMaxLineDown;
            windowAPI.minPatch.handle    @minMaxPatchDown;
            windowAPI.maxPatch.handle    @minMaxPatchDown;
            windowAPI.bigPatch.handle    @bigPatchDown
            };

        for k = 1 : size(buttonDownTable,1)
            h = buttonDownTable{k,1};
            callback = buttonDownTable{k,2};
            set(h, 'ButtonDownFcn', callback);
        end

        draggableObjList = [buttonDownTable{1:end-1,1}];
        cbk_id_cell = initCursorChangeOverDraggableObjs(hFig, draggableObjList);

        %====================================
        function minMaxLineDown(src,varargin)

            if src == windowAPI.maxLine.handle
                isMaxLine = true;
            else
                isMaxLine = false;
            end
            
            idButtonMotion = iptaddcallback(hFig, 'WindowButtonMotionFcn', ...
                                            @minMaxLineMove);
            idButtonUp = iptaddcallback(hFig, 'WindowButtonUpFcn', ...
                @minMaxLineUp);
            
            % Disable pointer manager.
            iptPointerManager(hFig, 'disable');

            %==============================
            function minMaxLineUp(varargin)

                acceptChanges(idButtonMotion, idButtonUp);
            end

            %====================================
            function minMaxLineMove(src,varargin)

                xpos = getCurrentPoint(hHistAx);
                if isMaxLine
                    newMax = xpos;
                    newMin = windowAPI.minLine.get();
                else
                    newMin = xpos;
                    newMax = windowAPI.maxLine.get();
                end
                newClim = validateClim([newMin newMax]);
                if isequal(newClim(1), xpos) || isequal(newClim(2), xpos)
                    updateAll(newClim);
                end
            end
        end %lineButtonDown

        %=================================
        function centerPatchDown(varargin)

            idButtonMotion = iptaddcallback(hFig, 'WindowButtonMotionFcn', ...
                                            @centerPatchMove);
            idButtonUp = iptaddcallback(hFig, 'WindowButtonUpFcn', @centerPatchUp);

            % Disable pointer manager.
            iptPointerManager(hFig, 'disable');

            startX = getCurrentPoint(hHistAx);
            oldCenterX = windowAPI.centerLine.get();

            %===============================
            function centerPatchUp(varargin)
                
                acceptChanges(idButtonMotion, idButtonUp);
            end

            %=================================
            function centerPatchMove(varargin)

                newX = getCurrentPoint(hHistAx);
                delta = newX - startX;

                % Set the window endpoints.
                centerX = oldCenterX + delta;
                minX = windowAPI.minLine.get();
                maxX = windowAPI.maxLine.get();
                width = maxX - minX;
                [newMin, newMax] = computeClim(width, centerX);
                newClim = validateClim([newMin newMax]);
                updateAll(newClim);
            end
        end %centerPatchDown

        %======================================
        function minMaxPatchDown(src, varargin)

            if isequal(src, windowAPI.minPatch.handle)
                srcLine = windowAPI.minLine;
                minPatchMoved = true;
            else
                srcLine = windowAPI.maxLine;
                minPatchMoved = false;
            end

            startX = getCurrentPoint(hHistAx);
            oldX = srcLine.get();
            
            idButtonMotion = iptaddcallback(hFig, 'WindowButtonMotionFcn', ...
                                            @minMaxPatchMove);
            idButtonUp = iptaddcallback(hFig, 'WindowButtonUpFcn', ...
                @minMaxPatchUp);

            % Disable pointer manager.
            iptPointerManager(hFig, 'disable');

            %===============================
            function minMaxPatchUp(varargin)

                acceptChanges(idButtonMotion, idButtonUp);
            end

            %======================================
            function minMaxPatchMove(src, varargin)

                newX = getCurrentPoint(hHistAx);
                delta = newX - startX;

                % Set the window endpoints.
                if minPatchMoved
                    minX = oldX + delta;
                    maxX = windowAPI.maxLine.get();
                else
                    maxX = oldX + delta;
                    minX = windowAPI.minLine.get();
                end
                newClim = validateClim([minX maxX]);
                updateAll(newClim);
            end
        end %minMaxPatchDown

        %==============================
        function bigPatchDown(varargin)

            idButtonMotion = iptaddcallback(hFig, 'windowButtonMotionFcn', ...
                                            @bigPatchMove);
            idButtonUp = iptaddcallback(hFig, 'WindowButtonUpFcn', @bigPatchUp);

            % Disable pointer manager.
            iptPointerManager(hFig, 'disable');

            startX = get(hHistAx, 'CurrentPoint');
            oldMinX = windowAPI.minLine.get();
            oldMaxX = windowAPI.maxLine.get();

            %============================
            function bigPatchUp(varargin)
                
                acceptChanges(idButtonMotion, idButtonUp);
            end

            %===========================
            function bigPatchMove(varargin)

                newX = getCurrentPoint(hHistAx);
                delta = newX(1) - startX(1);

                % Set the window endpoints.
                newMin = oldMinX + delta;
                newMax = oldMaxX + delta;

                % Don't let window shrink when dragging the window patch.
                origWidth = getWidthOfWindow;
                histRange = histStruct.histRange;
                
                if newMin < histRange(1)
                    newMin = histRange(1);
                    newMax = newMin + origWidth;
                end

                if newMax > histRange(2)
                    newMax = histRange(2);
                    newMin = newMax - origWidth;
                end
                newClim = validateClim([newMin newMax]);
                updateAll(newClim);
            end
        end %bigPatchDown
    
        %=================================================
        function acceptChanges(idButtonMotion, idButtonUp)
            
           iptremovecallback(hFig, 'WindowButtonMotionFcn', idButtonMotion);
           iptremovecallback(hFig, 'WindowButtonUpFcn', idButtonUp);
           
           % Enable the figure's pointer manager.
           iptPointerManager(hFig, 'enable');
           
           updateAllAndSaveInHistory(getClim);
           
        end
        
        %================================
        function width = getWidthOfWindow
            width = editBoxAPI.widthEdit.get();
        end

    end % setUpCallbacksOnDraggableWindow

    %=====================================
    function setUpCallbacksOnWindowWidgets

        callbackTable = {
            editBoxAPI.centerEdit  @actOnCenterChange;
            editBoxAPI.widthEdit   @actOnWidthChange;
            editBoxAPI.maxEdit     @actOnMinMaxChange;
            editBoxAPI.minEdit     @actOnMinMaxChange;
            };
        
        for m = 1 : size(callbackTable,1)
            h = callbackTable{m,1}.handle;
            callback = callbackTable{m,2};
            set(h, 'Callback', callback);
        end

        eyedropperAPI = clipPanelAPI.eyedropperAPI;
        droppers = [eyedropperAPI.minDropper.handle ...
                    eyedropperAPI.maxDropper.handle]; 
        set(droppers, 'callback', @eyedropper);

        %===================================
        function actOnMinMaxChange(varargin)

            areEditBoxStringsValid = checkEditBoxStrings;
            if areEditBoxStringsValid
                newMax = editBoxAPI.maxEdit.get();
                newMin = editBoxAPI.minEdit.get();
                [newClim] = validateClim([newMin newMax]);
                updateAllAndSaveInHistory(newClim);
            else
                resetEditValues;
                return;
            end
        end

        %==================================
        function actOnWidthChange(varargin)

            areEditBoxStringsValid = checkEditBoxStrings;
            if areEditBoxStringsValid
                centerValue = editBoxAPI.centerEdit.get();
                widthValue = editBoxAPI.widthEdit.get();

                [newMin newMax] = computeClim(widthValue, centerValue);
                newClim = validateClim([newMin newMax]); 
                updateAllAndSaveInHistory(newClim);
            else
                resetEditValues;
                return;
            end
        end

        %===================================
        function actOnCenterChange(varargin)

            areEditBoxStringsValid = checkEditBoxStrings;
            if areEditBoxStringsValid
                centerValue = editBoxAPI.centerEdit.get();
                widthValue = editBoxAPI.widthEdit.get();
                [newMin newMax] = computeClim(widthValue, centerValue);
                XLim = get(hHistAx,'XLim');

                % React to a center change that makes the newMin or 
                % newMax go outside of the XLim, but keep the center 
                % that the user requested.
                if ((newMin < XLim(1)) && (newMax > XLim(2)))
                    newMin = XLim(1);
                    newMax = XLim(2);
                elseif (newMin < XLim(1))
                    newMin = XLim(1);
                    newMax = newMin + 2 * (centerValue - newMin);
                elseif (newMax > XLim(2))
                    newMax = XLim(2);
                    newMin = newMax - 2 * (newMax - centerValue);
                end
                newClim = validateClim([newMin newMax]);
                updateAllAndSaveInHistory(newClim);
            else
                resetEditValues;
                return;
            end
        end

        %=======================
        function resetEditValues

            Clim = getClim;
            for k = 1 : size(callbackTable,1)
                callbackTable{k,1}.set(Clim);
            end
        end

        %=================================
        function eyedropper(src, varargin)

            if isequal(src, eyedropperAPI.minDropper.handle)
                editBox = editBoxAPI.minEdit;
                dropper = eyedropperAPI.minDropper;
            else
                editBox = editBoxAPI.maxEdit;
                dropper = eyedropperAPI.maxDropper;
            end

            % Prevent uicontrols from issuing callbacks before dropper is done.
            parent = ancestor(editBox.handle, 'uiflowcontainer', 'toplevel');
            children = findall(parent, 'Type', 'uicontrol');
            origEnable = get(children, 'Enable');
            set(children, 'Enable', 'off');

            % W/L mouse action sometimes conflicts afterward.  Turn it off briefly.
            origBDF = get(imageHandle, 'ButtonDownFcn');
            set(imageHandle, 'ButtonDownFcn', '');

            % Change the pointer to an eyedropper over the image.
            origPointerBehavior = iptGetPointerBehavior(imageHandle);
            enterFcn = @(f,cp) set(f, 'Pointer', 'custom', ...
                                      'PointerShapeCData', ...
                                      getEyedropperPointer(dropper.get()), ...
                                      'PointerShapeHotSpot', [16 1]);
            iptSetPointerBehavior(imageHandle, enterFcn);

            % Change the status text.
            origMsg = get(hStatusLabel, 'string');
            newMsg1 = 'Click on a pixel in the image to set the ';
            newMsg2 = sprintf('window''s %s value.', dropper.get());
            set(hStatusLabel, 'string', sprintf('%s%s',newMsg1,newMsg2));
            set(hStatusLabel, 'Enable', 'on');
            % Take care to undo all of these actions if the 
            % adjustment tool closes.
            origCloseRequestFcn = get(hFig, 'CloseRequestFcn');
            set(hFig, 'CloseRequestFcn', @closeDuringEyedropper)

            value = graysampler(imageHandle);

            % Set the edit text box.
            if (~isempty(value))
                editBox.set([value value]);
                areValid = checkEditBoxStrings;

                if areValid
                    newClim = [editBoxAPI.minEdit.get(), ...
                               editBoxAPI.maxEdit.get()];
                    newClim = validateClim(newClim);
                    updateAllAndSaveInHistory(newClim);
                else
                    resetEditValues;
                end
            end

            undoEyedropperChanges;
            
            %=====================================================
            function PointerShapeCData = getEyedropperPointer(tag)

                iconRoot = ipticondir;
                if strcmp(tag,'minimum')
                    cursor_filename = fullfile(iconRoot, ...
                                               'cursor_eyedropper_black.png');
                else
                    cursor_filename = fullfile(iconRoot, ...
                                               'cursor_eyedropper_white.png');
                end

                cdata = makeToolbarIconFromPNG(cursor_filename);
                PointerShapeCData = cdata(:,:,1)+1;
            end

            %=============================
            function undoEyedropperChanges

                % Change the pointer back.
                if ishandle(imageHandle)
                    iptSetPointerBehavior(imageHandle, origPointerBehavior);
                    
                    % Force pointer manager update.
                    iptPointerManager(ancestor(imageHandle, 'figure'));
                end

                % Change the message back.
                if ishandle(hStatusLabel)
                    set(hStatusLabel, 'string', origMsg);
                end

                % Turn the W/L mouse action back on if necessary.
                if ishandle(imageHandle)
                    set(imageHandle, 'ButtonDownFcn', origBDF);
                end

                % Reenable other uicontrols.
                for p = 1:numel(origEnable)
                    if ishandle(children(p))
                        set(children(p), 'Enable', origEnable{p});
                    end
                end
            end

            %=======================================
            function closeDuringEyedropper(varargin)

                undoEyedropperChanges;
                if ((~isempty(origCloseRequestFcn)) && ...
                        (~isequal(origCloseRequestFcn, 'closereq')))
                    feval(origCloseRequestFcn);
                end

                if (ishandle(hFig))
                    delete(hFig)
                end
            end

        end %eyedropper
        
        %======================================
        function areValid = checkEditBoxStrings

            centerValue = editBoxAPI.centerEdit.get();
            maxValue    = editBoxAPI.maxEdit.get();
            minValue    = editBoxAPI.minEdit.get();
            widthValue  = editBoxAPI.widthEdit.get();

            areValid = true;

            % Validate data.
            % - If invalid: display dialog, reset to last good value, stop.
            % - If valid: go to other callback processor.
            isValueEmpty = any([isempty(minValue), isempty(maxValue),...
                isempty(widthValue), isempty(centerValue)]);

            isValueString = any([ischar(minValue), ischar(maxValue),...
                ischar(widthValue), ischar(centerValue)]);

            isValueNonScalar = (numel(minValue) + numel(maxValue) +...
                numel(widthValue) + numel(centerValue) ~= 4);

            if (isValueEmpty || isValueString || isValueNonScalar)

                areValid = false;
                errordlg({'This window value must be numeric and scalar.'}, ...
                    'Invalid window value', ...
                    'modal')

            elseif (minValue >= maxValue)

                areValid = false;
                errordlg({'The minimum value must be less than maximum value.'}, ...
                    'Invalid window value', ...
                    'modal')

            elseif (((widthValue < 1) && (~isDoubleOrSingleData)) || ...
                    (widthValue <= 0))

                areValid = false;
                errordlg({'The window width must be greater than zero.'}, ...
                    'Invalid window value', ...
                    'modal')

            elseif ((floor(centerValue * 2) ~= centerValue * 2) && (~isDoubleOrSingleData))

                areValid = false;
                errordlg({'The window center value must be an integer.'}, ...
                    'Invalid window value', ...
                    'modal')
            end
        end %validateEditBoxStrings

    end %setUpCallbacksOnWindowWidgets

    %=========================================================
    function [minPixel, maxPixel] = computeClim(width, center)
        %FINDWINDOWENDPOINTS   Process window and level values.
        
        minPixel = (center - width/2);
        maxPixel = minPixel + width;
    end

    %===================================
    function setUpCallbacksOnAutoScaling
    
        callbackTable = {
            scalePanelAPI.elimRadioBtn       @changeScaleDisplay;
            scalePanelAPI.matchDataRangeBtn  @changeScaleDisplay;
            scalePanelAPI.scaleDisplayBtn    @autoScaleApply
            scalePanelAPI.percentEdit        @autoScaleApply;
        };
        
        for k = 1 : size(callbackTable,1)
            h = callbackTable{k,1}.handle;
            callback = callbackTable{k,2};
            set(h,'Callback', callback);
        end
        
        set(scalePanelAPI.percentEdit.handle, ...
            'ButtonDownFcn', @changeScaleDisplay, ...
            'KeyPressFcn', @changeScaleDisplay);

        % make matchDataRangeBtn selected by default.
        scalePanelAPI.matchDataRangeBtn.set(true);
        scalePanelAPI.elimRadioBtn.set(false);
        
        %========================================
        function changeScaleDisplay(src, varargin)

            if isequal(src, scalePanelAPI.matchDataRangeBtn.handle)
                scalePanelAPI.matchDataRangeBtn.set(true);
                scalePanelAPI.elimRadioBtn.set(false);
            else
                scalePanelAPI.matchDataRangeBtn.set(false);
                scalePanelAPI.elimRadioBtn.set(true);
            end
        end

        %================================
        function autoScaleApply(varargin)

            % Verify the percent and use it if box is checked.
            outlierPct = scalePanelAPI.percentEdit.get();

            matchDataRange = ...
                isequal(scalePanelAPI.matchDataRangeBtn.get(), true);

            CData = get(imageHandle, 'CData');
            minCData = min(CData(:));
            maxCData = max(CData(:));

            if matchDataRange

                localNewClim = [double(minCData) double(maxCData)];
                
            else
                % eliminate Outliers. 
                if isempty(outlierPct) || outlierPct > 100 || outlierPct < 0
                    errordlg({'Percentage must be a number between 0 and 100.'}, ...
                        'Invalid Percentage', ...
                        'modal')
                    scalePanelAPI.percentEdit.set('2');
                    return;
                end

                outlierPct = outlierPct / 100;

                % Double image data not in default range must be scaled and
                % shifted to the range [0,1] for STRETCHLIM to do 
                % the right thing.
                doubleImageOutsideDefaultRange = isDoubleOrSingleData && ...
                    (minCData < 0 || maxCData > 1);

                if doubleImageOutsideDefaultRange
                    % Keep track of old CData range for reconversion.
                     CData = mat2gray(CData);
                end

                localNewClim = stretchlim(CData, outlierPct / 2);

                if isequal(localNewClim, [0;1])
                    if outlierPct > 0.02
                        errordlg({'The specified percentage is too great.', ...
                            'Choose a smaller percentage.'}, ...
                            'Percentage Too Large', ...
                            'modal')
                        return;
                    elseif outlierPct ~= 0
                        errordlg({'This image contains too few grayscale values to eliminate outliers.',...
                            'Use the match data range option.'},...
                            'Cannot Eliminate Outliers',...
                            'modal')
                         return;
                    end
                end
                   
                % Scale the Clim from STRETCHLIM's [0,1] to match the range
                % of the data.
                if ~isDoubleOrSingleData
                    imgClass = class(CData);
                    localNewClim = double(intmax(imgClass)) * localNewClim;
                elseif doubleImageOutsideDefaultRange
                    localNewClim = localNewClim * (maxCData - minCData);
                    localNewClim = localNewClim + minCData;
                end
            end

            newClim = validateClim(localNewClim);
            updateAllAndSaveInHistory(newClim);

        end % autoScaleApply
    end %setUpCallbacksOnAutoScaling

    %====================================
    function newClim = validateClim(clim)

        % Prevent new endpoints from exceeding the min and max of the
        % histogram range, which is a little less than the xlim endpoints.
        % Don't want to get to the actual endpoints because there is a
        % problem with the painters renderer and patchs at the edge
        % (g298973).  histStruct is a variable calculated in the beginning
        % of createHistogramPalette.
        histRange = histStruct.histRange;
        newMin = max(clim(1), histRange(1));
        newMax = min(clim(2), histRange(2));
            
        if ~isDoubleOrSingleData
            % If the image has an integer datatype, don't allow the new endpoints
            % to exceed the min or max of that datatype.  For example, We don't
            % want to allow this because it wouldn't make sense to set the clim
            % of a uint8 image beyond 255 or less than 0.
            minOfDataType = double(intmin(getClassType(imgModel)));
            maxOfDataType = double(intmax(getClassType(imgModel)));
            newMin = max(newMin, minOfDataType);
            newMax = min(newMax, maxOfDataType);
        end
        
        % Keep min < max
        if ( ((newMax - 1) < newMin) && ~isDoubleOrSingleData )

            % Stop at limiting value.
            Clim = getClim;
            newMin = Clim(1);
            newMax = Clim(2);

            %Made this less than or equal to as a possible workaround to g226780
        elseif ( (newMax <= newMin) && isDoubleOrSingleData )

            % Stop at limiting value.
            Clim = getClim;
            newMin = Clim(1);
            newMax = Clim(2);
        end

        newClim = [newMin newMax];
    end


    %===================================================================
    function [hStatusLabel panelWidth] = createStatusPanel(parent)

        hStatusPanel = uipanel('Parent', parent, ...
            'Units', 'pixels', ...
            'Tag', 'StatusPanel',...
            'BorderType', 'none');

        defaultMessage = sprintf('%s%s', ...
            'Adjust the histogram above, or click ', ...
            'and drag the mouse over the image.');

        hStatusLabel = uicontrol('parent', hStatusPanel, ...
            'units', 'pixels', ...
            'style', 'text', ...
            'HorizontalAlignment', 'left', ...
            'string', defaultMessage);

        labelExtent = get(hStatusLabel, 'extent');
        panelWidth = labelExtent(3);
        panelHeight = labelExtent(4);
        set(hStatusLabel, 'Position', [1 1 panelWidth panelHeight]);
        set(hStatusPanel, 'HeightLimits', ...
                          [panelHeight panelHeight]);
                          
    end % createStatusPanel

    %======================
    function updateEditMenu
    
        % enable the undo menus when the clim gets its first change
        if currentHistoryIndex == 2
            set([undoMenu, undoAllMenu], 'Enable', 'on');
        elseif currentHistoryIndex == 1
            set([undoMenu, undoAllMenu], 'Enable', 'off');
        end

        % enable the redo menu when the length of the history is greater
        % than the current index
        historyLength = size(climHistory, 1);
        if historyLength > currentHistoryIndex
            set(redoMenu, 'Enable', 'on');
        elseif historyLength == currentHistoryIndex
            set(redoMenu, 'Enable', 'off');
        end
    end % updateEditMenu

    %===============================
    function undoLastChange(obj,evt)
        currentHistoryIndex = max(1,  currentHistoryIndex - 1);
        updateAll(climHistory(currentHistoryIndex,:));
        updateEditMenu
    end

    %=============================
    function redoLastUndo(obj,evt)
        historyLength = size(climHistory, 1);
        currentHistoryIndex = min(historyLength, currentHistoryIndex + 1);
        updateAll(climHistory(currentHistoryIndex,:));
        updateEditMenu
    end

    %===============================
    function undoAllChanges(obj,evt)
        currentHistoryIndex = 1;
        updateAll(climHistory(currentHistoryIndex,:));
        updateEditMenu
    end

    %==========================================
    function updateAllAndSaveInHistory(newClim)
        % get the length of entries in the history
        historyLength = size(climHistory,1);

        % increment current index by one to indicate the new entry's
        % position.
        currentHistoryIndex = currentHistoryIndex + 1;

        % if the length of entries in the history is longer that the
        % current index we discard all entries after the current index.
        if historyLength > currentHistoryIndex
            climHistory(currentHistoryIndex,:) = [];
        end
        climHistory(currentHistoryIndex,:) = [newClim(1), newClim(2)];

        updateAll(newClim);
        updateEditMenu;
    end

    %==========================
    function updateAll(newClim)

        % Update edit boxes with new values.
        updateEditBoxes(newClim);

        % Update patch display.
        updateHistogram(newClim);

        % we don't want the clim event handler executed to prevent
        % duplicate entries in the history.
        blockEventHandler = true;

        % Update image Clim.
        updateImage(hImageAx, newClim);

        blockEventHandler = false;
    end

    %===============================
    function updateEditBoxes(newClim)
    
        names = fieldnames(editBoxAPI);
        for k = 1 : length(names)
            editBoxAPI.(names{k}).set(newClim);
        end
    end 

    %================================
    function updateHistogram(newClim)

        names = fieldnames(windowAPI);
        for k = 1 : length(names)
            windowAPI.(names{k}).set(newClim);
        end
    end % updateHistogram

    %===================================
    function updateImage(hImageAx, clim)

        if clim(1) >= clim(2)
            eid = sprintf('Images:%s:internalError',mfilename);
            msg = 'Internal error - clim(1) is >= clim(2).';
            error(eid,'%s',msg);
        end
        set(hImageAx, 'clim', clim);
    end

    %======================
    function clim = getClim
        clim = get(hImageAx,'Clim');
    end

end % createHistogramPalette


%==========================================================================
function imgModel = validateImage(hIm)

imgModel = getimagemodel(hIm);
if ~strcmp(getImageType(imgModel),'intensity')
  eid = sprintf('Images:%s:unsupportedImageType',mfilename);
  error(eid, 'Only intensity images are supported.');
end

cdata = get(hIm,'cdata');
if isempty(cdata)
    eid = sprintf('Images:%s:invalidImage',mfilename);
    error(eid, '%s does not support empty images.',upper(mfilename));
end

end

%==========================================================================
function cbk_id_cell = initCursorChangeOverDraggableObjs(client_fig, drag_objs)
% initCursorChangeOverDraggableObjs

% initialize variables for function scope
num_of_drag_objs    = numel(drag_objs);

enterFcn = @(f,cp) setptr(f, 'lrdrag');
iptSetPointerBehavior(drag_objs, enterFcn);

% Add callback to turn on flag indicating that dragging has stopped.
stop_drag_cbk_id = iptaddcallback(client_fig, ...
    'WindowButtonUpFcn', @stopDrag);

obj_btndwn_fcn_ids = zeros(1, num_of_drag_objs);

% Add callback to turn on flag indicating that dragging has started
for n = 1 : num_of_drag_objs
    obj_btndwn_fcn_ids(n) = iptaddcallback(drag_objs(n), ...
        'ButtonDownFcn', @startDrag);
end

cbk_id_cell = {client_fig, 'WindowButtonUpFcn', stop_drag_cbk_id;...
    drag_objs,  'ButtonDownFcn', obj_btndwn_fcn_ids};


    %==========================
    function startDrag(obj,evt)
        % Disable the pointer manager while dragging.
        iptPointerManager(client_fig, 'disable');
    end

    %========================
    function stopDrag(ob,evt)
        % Enable the pointer manager.
        iptPointerManager(client_fig, 'enable');
    end

end % initCursorChangeOverDraggableObjs


%==========================================================================
function deleteCursorChangeOverDraggableObjs(cbk_id)

rows = size(cbk_id);
for n = 1 : rows
    id_length = length(cbk_id{n,1});
    for m = 1 : id_length
        iptremovecallback(cbk_id{n,1}(m), cbk_id{n,2}, cbk_id{n,3}(m));
    end
end
end % deleteCursorChangeOverDraggableObjs


%===========================================================================
function wdlg = displayWarningDialog(curClim, imDataLim)

str1 = 'You cannot use the Adjust Contrast tool with the existing';
str2 = 'display range because the limits of the display range [%s %s]';
str3 = 'do not fall inside the limits of the data range [%s %s].';
tmpstr1 = sprintf('%s %s %s', str1, str2, str3);

str4 = '\nClick OK to adjust the display range so that you can use the tool.';

formatValue = @(v) sprintf('%0.0f', v);
str{1} = sprintf(tmpstr1, ...
    formatValue(curClim(1)), formatValue(curClim(2)), ...
    formatValue(imDataLim(1)), formatValue(imDataLim(2)));

str{2} = sprintf(str4);
wdlg = questdlg(str, ...
    'Invalid Display Range', 'OK', ...
    'Cancel', 'OK');
end

    %==========================================================
    function badValue = isDisplayRangeOutsideDataRange(him,hax)

        % Checking to see if the display range is outside the image's data range.
        clim = get(hax,'Clim');
        histStruct = getHistogramData(him);
        histRange = histStruct.histRange;
        badValue = false;

        if clim(1) < histRange(1) || clim(2) > histRange(2)
            badValue = true;
        end
    end

