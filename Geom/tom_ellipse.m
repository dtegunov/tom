function [out,mjr,mnr]=tom_ellipse(varargin)
%ELLIPSE Ellipse Grapics Object Using Line or Patch.
% He = ELLIPSE(Axes,Position,Angle) draws an ellipse using a line object on
% the current axes where:
%
% Axes = [Major Minor] is a vector containing the semi-minor axis length
% and semi-major axis length, where the semi-minor axis length is the
% distance between the ellipse center and the closest point on the ellipse
% and the semi-major axis length is the distance between the ellipse center
% and the farthest point on the ellipse.
%
% Position = [Xcenter Ycenter] is a vector containing the ellipse center,
%
% Angle is the angle in degrees of the Major axis with respect to the X axis,
%
% He is the handle of the graphics object created.
%
% He = ELLIPSE(Axes,Position,Angle,'PName',PValue,...) in addition sets
% ellipse properties using the provided property name-values pairs.
%
% He = ELLIPSE(Axes,Position,Angle,S) in addition sets ellipse properties
% using the provided property name-values pairs contained in the structure
% S, whose field names are property names and contents are associated
% property values.
%
% He = ELLIPSE('PName',PValue,...) or He = ELLIPSE(S) creates an ellipse
% using only property name-value pairs.
%
% [He,MJR,MNR] = ELLIPSE(...) in addition returns handles to minor and
% major axis line objects respectively.
%
% ELLIPSE(He,'PName',PValue,...) or ELLIPSE(He,S) sets the properties of
% the ellipse having handle He using the provided property name-value
% pairs.
%
% ELLIPSE(He,'PName') returns the property value associated with the
% property name 'PName' for the ellipse having handle He.
%
% S = ELLIPSE(He) returns a structure S whose field names are all the
% property names and contents are all the associated property values.
%
% Properties:
% NAME         VALUE {default}
% Axes         [Major Minor], same as Axes input variable
% Position     [Xcenter Ycenter], same as Position input variable
% Angle        [degrees], same as Angle input variable
% Type         {'line'} or 'patch', read only after creation
% Ellipse      {'on'} or 'off' visibility of the ellipse
% Major        'on' or {'off'} visibility of a line marking the Major axis
% Minor        'on' or {'off'} visibility of a line marking the Minor axis
% Handles      [MJR;MNR], read only two element vector containing the
%              handles of the major and minor axis lines respectively
%
% Properties of the line or patch ellipse object created can be manipulated
% by using set and get with the returned handle He.
% Properties of the major and minor axis line objects can be manipulated by
% using set and get with the line object handles MJR and MNR.
%
% 'axis equal' must be applied to correct axis distortion.

% D.C. Hanselman, University of Maine, Orono, ME  04469-5708
% Mastering MATLAB 7
% www.eece.maine.edu/mm
% revised 2006-03-16, 2006-06-27, 2006-06-30

% parse inputs
if nargin==0
   error('At Least One Input Argument Required.')
elseif numel(varargin{1})==1 && ishandle(varargin{1})     % ELLIPSE(He,...)
   he=varargin{1};
   Sdef=getappdata(he,'Ellipse');
   if nargin==1                                     % ellipse(He), return S
      out=Sdef;
      return
   elseif nargin==2
      if isstruct(varargin{2})                         % ellipse(He,S), set
         S=local_get_S(Sdef,varargin{2});
      elseif ischar(varargin{2})       % ellipse(He,'PName'), return PValue
         pname=local_get_pname(varargin{2});
         out=Sdef.(pname);
         return
      else
         error('Unknown Second Argument.')
      end
   else                               % ellipse(He,'Pname',PValue,...), set
      S=local_get_S(Sdef,varargin{2:end});
   end
   %            set properties set properties set properties set properties
   local_check_validity(S)
   if ~isequal(S.Handles,Sdef.Handles)
      error('Handles Property is Read Only.')
   end
   if ~isequal(S.Type,Sdef.Type)
      error('Type Property is Read Only.')
   end
   theta=linspace(0,2*pi,500);
   Xdata=max(S.Axes)*cos(theta);
   Ydata=min(S.Axes)*sin(theta);
   Ydata(end)=0;
   Xdata(end)=max(S.Axes);
   d=S.Angle*pi/180;
   x=cos(d)*Xdata-sin(d)*Ydata + S.Position(1);
   y=sin(d)*Xdata+cos(d)*Ydata + S.Position(2);   
   set(he,'Xdata',x,'Ydata',y,'Visible',S.Ellipse);               % ellipse
   
   x=min(S.Axes)*cos(d+pi/2)*[-1 1] + S.Position(1);
   y=min(S.Axes)*sin(d+pi/2)*[-1 1] + S.Position(2);
   set(S.Handles(2),'Xdata',x,'Ydata',y,'Visible',S.Minor);     % axis line
   
   x=max(S.Axes)*cos(d)*[-1 1] + S.Position(1);
   y=max(S.Axes)*sin(d)*[-1 1] + S.Position(2);
   set(S.Handles(1),'Xdata',x,'Ydata',y,'Visible',S.Major);     % axis line
   
   setappdata(he,'Ellipse',S)

else                                                     % ellipse creation

   S=struct('Axes',[1 1],'Position',[0 0],'Angle',0,'Type','line',...
      'Ellipse','on','Major','off','Minor','off','Handles',[nan nan]);
   pidx=find(cellfun(@(x) ischar(x)||isstruct(x),varargin),1);
   % the above line requires version 7.1 or later
   % replace the above line with this for loop for earlier versions
   %---------------------------------------------------
%    pidx=[];
%    for k=1:length(varargin)
%       if ischar(varargin{k}) || isstruct(varargin{k})
%          pidx=k;
%          break
%       end
%    end
   %---------------------------------------------------
   if isempty(pidx)
      pidx=nargin+1;
   end
   if pidx>1
      S.Axes=varargin{1};
   end
   if pidx>2
      S.Position=varargin{2};
   end
   if pidx>3
      S.Angle=varargin{3};
   end
   if pidx<nargin+1
      if ischar(varargin{pidx})               % ellipse(...,'PName',PValue)
         S=local_get_S(S,varargin{pidx:end});
      elseif isstruct(varargin{pidx})                      % ellipse(...,S)
         S=local_get_S(S,varargin{pidx});
      else
         error('Unknown Input Syntax.')
      end
   end
   theta=linspace(0,2*pi,500);
   Xdata=max(S.Axes)*cos(theta);
   Ydata=min(S.Axes)*sin(theta);
   Ydata(end)=0;
   Xdata(end)=max(S.Axes);
   d=S.Angle*pi/180;
   x=cos(d)*Xdata-sin(d)*Ydata + S.Position(1);
   y=sin(d)*Xdata+cos(d)*Ydata + S.Position(2);
   if strncmpi(S.Type,'l',1)                                      % ellipse
      he=line(x,y,'Color','r','Visible',S.Ellipse);
   else
      he=patch('Xdata',x,'Ydata',y,...
         'EdgeColor','r','FaceColor','none','Visible',S.Ellipse);
   end
   x=min(S.Axes)*cos(d+pi/2)*[-1 1] + S.Position(1);
   y=min(S.Axes)*sin(d+pi/2)*[-1 1] + S.Position(2);
   mnr=line('Xdata',x,'Ydata',y,...
                     'Color','r','Visible',S.Minor);      % Minor axis line
   
   x=max(S.Axes)*cos(d)*[-1 1] + S.Position(1);
   y=max(S.Axes)*sin(d)*[-1 1] + S.Position(2);
   mjr=line('Xdata',x,'Ydata',y,...
                     'Color','r','Visible',S.Major);      % Major axis line

   S.Handles=[mjr;mnr];

   setappdata(he,'Ellipse',S)
   if nargout  % only return handle if asked
      out=he;
   end  
end
%--------------------------------------------------------------------------
function S=local_get_S(S,varargin)
if nargin==2 && isstruct(varargin{1})                      % ellipse(...,S)
   Snew=varargin{1};
   Snames=fieldnames(Snew);
   for k=1:length(Snames)
      pname=local_get_pname(Snames{k});
      S.(pname)=Snew.(Snames{k});
   end
elseif rem(length(varargin),2)~=0 || ~iscellstr(varargin(1:2:end))
   error('PNames and PValues Must Appear in Pairs')
else                                      % ellipse(...,'PName',PValue,...)
   for k=1:2:length(varargin)-1
      pname=local_get_pname(varargin{k});
      S.(pname)=varargin{k+1};
   end
end
%--------------------------------------------------------------------------
function pname=local_get_pname(pname)
% get and check property name
PNames={'Axes','Position','Angle','Type',...
        'Ellipse','Major','Minor','Handles'};
idx=find(strncmpi(pname,PNames,length(pname)));
if isempty(idx)   % no matches
   error(['Unknown Property Name: ' pname])
elseif length(idx)>1 % more than one match
   error(['Property Name Not Unique: ' pname])
end
pname=PNames{idx};
%--------------------------------------------------------------------------
function local_check_validity(S)
if numel(S.Axes)~=2 || ~isnumeric(S.Axes)
   error('Axes Must be a Two Element Numeric Vector.')
end
if numel(S.Position)~=2 || ~isnumeric(S.Position)
   error('Position Must be a Two Element Numeric Vector.')
end
if numel(S.Angle)~=1 || ~isnumeric(S.Angle)
   error('Angle Must be a Numeric Scalar.')
end
if ~ischar(S.Type)
   error('Type Property Must be a String.')
elseif isempty(strncmpi(S.Type,{'line' 'patch'},1))
   error('Type Property Must be Either ''Line'' or ''Patch''.')
end
if ~ischar(S.Ellipse) || isempty(strncmpi(S.Ellipse,{'on' 'off'},2))
   error('Ellipse Property Must be Either ''On'' or ''Off''.')
end
if ~ischar(S.Major) || isempty(strncmpi(S.Major,{'on' 'off'},2))
   error('Major Property Must be Either ''On'' or ''Off''.')
end
if ~ischar(S.Minor) || isempty(strncmpi(S.Minor,{'on' 'off'},2))
   error('Minor Property Must be Either ''On'' or ''Off''.')
end
