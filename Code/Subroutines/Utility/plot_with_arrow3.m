function [plotHandle] = plot_with_arrow3(data, varargin)
% PLOT.WITH_ARROW3  Plot data with arrows along the curve showing the direction
%   This function is meant to extend the abilities of the PLOT3 command by
%   adding arrowheads to the curve.
%
%   H = plotWithArrows(data) - plot the specified data with the default
%   options. The function returns the handle to the plot
%   
%   H = plotWithArrows(data, LineType) - plot the specified data with
%   the specified line type (e.g. 'b*' or '--r'). See PLOT3 help for more
%   details
%
%   H = plotWithArrows(data, LineType, ...) - Allows you to specify
%   additional options for the plot. See the OPTIONS section below for
%   details.
%
%   OPTIONS:
%
%       'Color'         -   color string | rgb array - The color of the
%                           data and arrowheads. Default is 'b'
%
%       'LineWidth'     -   a number specifying the width of the plotted
%                           line. See PLOT3 help for more details. Default
%                           is 1.25
%
%       'DisplayName'   -   a string for legend
%
%       'HandleVisibility'  'on' | 'off' - visible in legend
%
%       'NumArrows'     -   The number of arrowheads to place on the data.
%                           Default is 4
%
%       'AlignArrows'   -   'center' | 'start' | 'end' - alignment of
%                           arrows along trajectory.
%                           Default is 'center'
%
%       'FlipDir'       -   False | True - Whether or not to flip the
%                           direction of the arrowheads (e.g. if the data
%                           is listed in reverse time). Default is false
%
%       'ArrowScale'    -   Scale factor for the arrowheads. Default is 1
%
%       'BaseColor'     -   Color on the base of arrowheads.
%                           Default scales Color value by 0.8
%
%       'NumFaces'      -   Number of faces per side of the arrow.
%                           Default is 16
%
%       'Offset'        -   Offsets arrows by a fraction of the line length
%   See also PLOT3
%
%   Author: Robert Lee
%   Version: November 14, 2024
%   Editted by: Jay Singh (November 20, 2024)
%

% TODO: add even time spacing

%% Set up Parsing
defLineWidth = 2;
defDisplayName = '';
defHandleVisibility = 'off';
defNumArr = 4;
defAlignArr = 'center';
defFlipDir = false;
defArrowScale = 1.5;
defLineType = '-';
defColor = gca().ColorOrder(gca().ColorOrderIndex,:);
defNumFaces = 16;
defOffset = 0;

p = inputParser;

validColor = @(x) (isnumeric(x) && length(x) == 3) || ischar(x) || isstring(x);
notNegNum = @(x) isnumeric(x) && x >= 0;

addRequired(p, 'data', @isnumeric);

addOptional(p, 'LineType', defLineType, @(x) ischar(x) || isstring(x));

addParameter(p, 'Color', defColor, validColor);
addParameter(p, 'LineWidth', defLineWidth, notNegNum);
addParameter(p, 'DisplayName', defDisplayName, @(x) ischar(x) || isstring(x));
addParameter(p, 'HandleVisibility', defHandleVisibility, @(x) ischar(x) || isstring(x));
addParameter(p, 'NumArrows', defNumArr, notNegNum);
addParameter(p, 'AlignArrows', defAlignArr, @(x) ischar(x) || isstring(x));
addParameter(p, 'FlipDir', defFlipDir, @islogical);
addParameter(p, 'ArrowScale', defArrowScale, notNegNum);
addParameter(p, 'NumFaces', defNumFaces, notNegNum);
addParameter(p, 'Offset', defOffset, @isnumeric);

%% Parse
% Check to see if the first optional input could be a LineType
if(~isempty(varargin) > 0 && ((ischar(varargin{1}) && length(varargin{1}) < 5) || (isstring(varargin{1}) && strlength(varargin{1}) < 5)) )
    % It is, proceed normally
    parse(p, data, varargin{:});
else
    % It isn't, throw in the default to avoid errors and parse the optional
    % inputs as param-value pairs
    parse(p, data, defLineType, varargin{:});
end
xdata = data(1,:);
ydata = data(2,:);
zdata = data(3,:);
color = p.Results.Color;
lineType = p.Results.LineType;
lineWidth = p.Results.LineWidth;
displayName = p.Results.DisplayName;
handleVisibility = p.Results.HandleVisibility;
numArrows = p.Results.NumArrows;
alignArrows = p.Results.AlignArrows;
flipDir = p.Results.FlipDir;
scale = p.Results.ArrowScale;
numFaces = p.Results.NumFaces;
offset = p.Results.Offset;

% Error checking
if(size(data,1) < 3)
    error('data does not have at least 3 rows');
end

if(numArrows > length(data))
    error('Number of arrows exceeds number of data points! Cannot create arrows...');
end


% Make the color match a color specified in LineType
% letters = ischarprop(lineType, 'alpha'); %logical array: whether or not each char is a letter
letters = isletter(lineType);
if(sum(letters) > 0)
    color = lineType(letters);
end

%% Compute arrow directions and locations
arrows = zeros(numArrows, numFaces+2, 3);
stepSize = floor(length(xdata)/numArrows);
if strcmpi(alignArrows,"center")
  stepStart = floor(length(xdata)/(2*numArrows))+1;
elseif strcmpi(alignArrows,"start")
  stepStart = 1;
elseif strcmpi(alignArrows,"end")
  stepStart = length(xdata) - (numArrows-1)*stepSize;
else
  error("Invalid arrow alignment specified");
end

% Range of x,y,z; use to choose size of arrowhead
xExtent = abs(max(xdata) - min(xdata));
yExtent = abs(max(ydata) - min(ydata));
zExtent = abs(max(zdata) - min(zdata));
avgExt = mean([xExtent, yExtent, zExtent]);

% Compute dimensions
l = 0.04*avgExt*scale;  % Length of arrowhead
w = 0.8*l;                  % Width of arrowhead
s = -0.6*l;             % Distance from base point to bottom (flat edge) of arrowhead
m = 0.2*l;              % Indent distance from bottom (flat edge) of arrowhead

offsetStep = round(offset*length(xdata));
for n = 1:numArrows
  ix = (n-1)*stepSize+stepStart+offsetStep;
  if(ix > length(xdata)); break; end
  
  loc = [xdata(ix), ydata(ix), zdata(ix)];
  
  if(ix < length(xdata))
    dir = [xdata(ix+1), ydata(ix+1), zdata(ix+1)] - loc;
  else
    dir = loc - [xdata(ix-1), ydata(ix-1), zdata(ix-1)];
  end
  
  % Normalize length of dir and flip it if desired
  dir = (-1)^flipDir * dir/norm(dir);
  
  % Generate perpendicular direction
  % source: https://math.stackexchange.com/a/4112622
  s_xz = sign((sign(dir(1))+0.5)*(sign(dir(3))+0.5));
  s_yz = sign((sign(dir(2))+0.5)*(sign(dir(3))+0.5));
  prp = [s_xz*dir(3), s_yz*dir(3), -s_xz*dir(1) - s_yz*dir(2)];
  prp = (w/2)*prp/norm(prp);

  arrows(n,1,:) = loc + (s+l)*dir;
  arrows(n,2,:) = loc + (s+m)*dir;
  % NumFaces points along base circle; use patch3() to fill these points later
  th = linspace(0,2*pi,numFaces);
  for i = 1:numFaces
    % Axis-angle rotation
    % source: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    arrows(n,i+2,:) = loc + s*dir ...
      + prp*cos(th(i)) + cross(dir,prp)*sin(th(i)) + dir*dot(dir,prp)*(1-cos(th(i)));
  end
end

%% Plotting

% are we holding?
held = ishold;

if(~ishold); hold on; end

figure(gcf);
plotHandle = plot3(xdata, ydata, zdata, lineType, ...
  'Color', color, ...
  'LineWidth', lineWidth, ...
  'HandleVisibility', handleVisibility, ...
  'DisplayName', displayName);
for n = 1:size(arrows,1)
  % draw bottom
  patch('Vertices', squeeze(arrows(n,:,:)), ...
    'Faces', [2,(1:numFaces)+2,2], ...
    'FaceColor', 0.6*color, ...
    'LineStyle', 'none', ...
    'HandleVisibility', 'off');
  % draw top
  patch('Vertices', squeeze(arrows(n,:,:)), ...
    'Faces', [1,(1:numFaces)+2,1], ...
    'FaceColor', color, ...
    'EdgeColor', color, ...
    'HandleVisibility', 'off');
end

if(held); hold on; else; hold off; end