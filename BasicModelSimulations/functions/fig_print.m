function [] = fig_print(FIG,dir_filename)

%% Input
% 1- figure handle
% 2- Full path including the name


figure(FIG)
% GCA = gca;
% 
% GCA.Box = 'off';
% GCA.TickDir = 'out';

% GCA.FontSize = ft_size;

FIG.Units = 'Inches';
pos = FIG.Position;
FIG.PaperPositionMode = 'Auto';
FIG.PaperUnits = 'Inches';
FIG.PaperSize = [pos(3),pos(4)];
print(gcf,'-painters',dir_filename,'-dpdf')
% print(gcf,'-painters',dir_filename,'-dsvg')
% print(gcf,'-painters',dir_filename,'-depsc')