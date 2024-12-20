function [] = savepdfs(varargin)

% Converts everything in the folder to pdf
close all
%set(0,'DefaultFigureVisible','off')
% Get graphics path
if nargin == 0
    fileloc = '';
else
    fileloc = varargin(1);
    fileloc = fileloc{1};
end
if ispc
    currentpath = pwd;
    gp = fullfile(erase(currentpath,'\Code\H1'), '/Solutions/Parts/H1/graphics',fileloc);
else
    currentpath = pwd;
    gp = fullfile(erase(currentpath,'/Code/H1'), '/Solutions/Parts/H1/graphics',fileloc);    
end
rmpath(genpath(fullfile(currentpath, 'Output Figures')))
addpath(fullfile(currentpath, 'Output Figures',fileloc));

% Get all files
files = dir(fullfile("Output Figures",fileloc,'*.fig'));
figures = struct2cell(files);
for i = 1:length(files)
    outputName = [erase(figures{1,i},'.fig'),'.pdf'];
    savefig = openfig(figures{1,i});
    saveloc = fullfile(gp,outputName);
    exportgraphics(savefig, saveloc,Resolution=200);
end
addpath(genpath(fullfile(currentpath, 'Output Figures')))
close all
%set(0,'DefaultFigureVisible','on')