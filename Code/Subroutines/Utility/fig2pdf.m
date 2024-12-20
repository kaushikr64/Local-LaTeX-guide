function [] = fig2pdf(file,GraphicsPath)
%SAVEFIG2PDF Take a guess
    savedFig = [file,'.fig'];
    openfig(savedFig);
    exportgraphics(gcf, fullfile(GraphicsPath,[file,'.pdf']),Resolution=150);
    close all
end

