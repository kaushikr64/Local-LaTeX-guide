function figure_settings(font_size)

if nargin < 1
    font_size = 'small';
end

%% font
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultConstantLineInterpreter', 'latex');
set(groot, 'defaultTextboxshapeInterpreter', 'latex');
set(groot, 'defaultTextarrowshapeInterpreter', 'latex');
set(0,'defaultfigurecolor',[1 1 1])
set(groot, 'defaultTextboxshapeFaceAlpha', 0.5)

% GUI font
% set(groot, 'defaultUicontrolFontName', 'Meiryo UI');
set(groot, 'defaultUicontrolFontName', 'Times');

% Axis font
set(groot, 'defaultAxesFontName', 'Times');

% Title, legend, constant line font
set(groot, 'defaultTextFontName', 'Times');
set(groot, 'defaultConstantLineFontName', 'Times');

% GUI fontsize
set(groot, 'defaultUicontrolFontSize', 9);

set(groot,'defaultTiledLayoutPadding', 'tight')
set(groot,'defaultTiledLayoutTileSpacing', 'tight')

% Axis fontsize
switch font_size
    case 'large'
        set(groot, 'defaultAxesFontSize', 25);
        set(groot, 'defaultLegendFontSize', 25);
        set(groot, 'defaultColorbarFontSize', 25);
        set(groot, 'defaultTextarrowshapeFontSize', 18);
        set(groot, 'defaultTextboxshapeFontSize', 18);
        set(groot, 'defaultConstantLineFontSize', 18);

        % Title, legend fontsize
        set(groot, 'defaultTextFontSize', 18);
    case 'medium'
        set(groot, 'defaultAxesFontSize', 20);
        set(groot, 'defaultLegendFontSize', 20);
        set(groot, 'defaultColorbarFontSize', 20);
        set(groot, 'defaultTextarrowshapeFontSize',20);
        set(groot, 'defaultTextboxshapeFontSize', 20);
        set(groot, 'defaultConstantLineFontSize', 20);

        % Title, legend fontsize
        set(groot, 'defaultTextFontSize', 20);

    case 'small'
        set(groot, 'defaultAxesFontSize', 15);
        set(groot, 'defaultLegendFontSize', 15);
        set(groot, 'defaultColorbarFontSize', 15);
        set(groot, 'defaultTextarrowshapeFontSize',15);
        set(groot, 'defaultTextboxshapeFontSize', 15);
        set(groot, 'defaultConstantLineFontSize', 15);

        % Title, legend fontsize
        set(groot, 'defaultTextFontSize', 15);
end


%% other
% axis line thickness
set(groot, 'DefaultAxesLineWidth', 1);
set(groot, 'DefaultLegendAutoUpdate', 'off');

% legend object line thickness
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 10);
set(groot, 'DefaultStairLineWidth', 1.5);
set(groot, 'DefaultStairMarkerSize', 10);

% constant line (xline, yline) color
set(groot, 'defaultConstantLineColor', 'k');

% default color map
% cmap = spring(128);
% set(groot, 'defaultFigureColormap', cmap);

%Plot colors
corder = [0,0,1;1,0,0;0.8,0.8,0;
    0.66015625,0.66015625,0.66015625;
    0,0,0;1,0.64453125,0;
    1,0,1;0,0.5,0.5;
    0,0,0.54296875;0,0.390625,0;
    0,1,1;0.59765625,0.1953125,0.796875];
set(0, 'defaultAxesColorOrder', corder);

% figure size
set(groot, 'defaultFigureUnits','pixels')
% set(groot, 'defaultFigurePosition',[100 100 650 400])
% set(groot, 'defaultFigurePosition',[100 100 400 400])
width = 20; %cm
hwratio = 0.65;
set(groot, 'defaultFigureUnits', 'centimeters')
set(groot, 'defaultFigurePosition', [3 3 width, hwratio*width])

% grid on
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot,'defaultAxesZGrid','on')

% legend settings
set(groot, 'defaultLegendBox', 'on')
% set(groot, 'defaultLegendItemTokenSize', [10, 10])

% GUI (figure) settings
% set(groot, 'defaultTextboxshapeEdgeColor', [1 1 1])
set(groot, 'defaultTextboxshapeLineStyle', 'none')

% axes clipping (looks better than the default '3dbox')
set(groot, 'defaultAxesClippingStyle', 'rectangle')


set(groot,'DefaultFigureWindowStyle','normal')
set(groot,'DefaultFigureColor',[1,1,1])

% Set OpenGL
opengl hardwarebasic


% use this command to get all graphics objects
% types = unique(get(findall(gcf, '-property', 'Type'), 'Type'));
