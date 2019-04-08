%% Date created 06.05.18 by M. Mohagheghi

% This script combines the two maps computed by rudimentary classification
% of rebound, disinh, blue regions to the rebound, disinh, and spontaneous
% activity by observing firing rate.

% clear
% close all

function [] = maps_combination_modcolormap_simple_contour(proc_data_fl, res_dir)

load(proc_data_fl)

% res_dir = './colormap-100Hz-5ms-exc/';

if exist(res_dir) ~= 7
    mkdir(res_dir)
end

reb_th = 0.01;
ch_th = reb_th*100;

sel_GS = G_SNr>0;
sel_GC = G_CX>0;

% G_SNr = G_SNr(sel_Gs);
% G_CX = G_CX(sel_Gs);
% reb_af_CaKO = reb_af_CaKO(sel_Gs,sel_Gs);
% reg_af_bef  = reg_af_bef(sel_Gs,sel_Gs);
% no_spk_af_mov = no_spk_af_mov(sel_Gs,sel_Gs);

reb_disinh_map = 3*ones(length(G_SNr),length(G_CX));
reb_disinh_map(reb_af_CaKO>=ch_th)=1;
reb_disinh_map(logical((reg_af_bef>=ch_th)-(reb_af_CaKO>=ch_th)))=2;
reb_disinh_map(no_spk_af_mov<ch_th)=0;

mycolormap = cat(1,[0,0,0],[1,0,0],[1,1,0],[0,0,1],[0,1,0],[1,0,1]);

%% Colormap
Aut = flipud(autumn(100));
Aut = jet(100);
Aut = parula(100);
% Spr = summer(100);

% num_cols = 3;
nan_c = [0.4,0.4,0.4];
% sel_hues = [60,120,180,240,300,360];
% Hue_vals = [];
% Val_vals = [];
% Sat_vals = [];
% num_steps = 100;
% 
% for nc_ind = 1:num_cols-1
%     Hue_vals = [Hue_vals;sel_hues(nc_ind)*ones(num_steps,1)];
% %     Hue_vals = [Hue_vals;linspace(nc_ind-1,nc_ind,num_steps)'];
%     Sat_vals = [Sat_vals;linspace(0.25,1,num_steps)'];
% %     Sat_vals = [Sat_vals;ones(num_steps,1)];
%     Val_vals = [Val_vals;ones(num_steps,1)];
% %     Val_vals = [Val_vals;linspace(0.75,1,num_steps)'];
% end
% 
% Hue_vals = Hue_vals/max(Hue_vals);
% 
% HSV = [Hue_vals,Sat_vals,Val_vals];
% RGB = hsv2rgb(HSV);
% RGB = [black;RGB];
RGB = [nan_c;Aut];
%%

figure;
colormap(mycolormap)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),reb_disinh_map(sel_GS,sel_GC))
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
% fig_print(gcf,[res_dir,'rebound-disinhibition-map'])

reb_combined = ones(length(G_SNr),length(G_CX));
reb_combined(isnan(contaminated_reb))=0;
reb_combined(contaminated_reb>=reb_th)=1;
reb_combined(contaminated_reb<reb_th)=2;

figure;
colormap(mycolormap)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),reb_combined(sel_GS,sel_GC))
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];

% nonz_inds = reb_combined == 0;
% comb_all_class = reb_disinh_map;
% comb_all_class = reb_disinh_map + reb_combined;

combine_all = zeros(length(G_SNr),length(G_CX));

% % Rebound
% combine_all(reb_disinh_map==1 & reb_combined==1) = 0;
% % combine_all = combine_all + reb_af_CaKO/max(reb_af_CaKO(:));
% combine_all(reb_disinh_map==0 & reb_combined==1) = 0;
% 
%Rebound + disinhibition
% combine_all(reb_disinh_map==2 & reb_combined==1) = 0;
cont_reb_lim = contaminated_reb;
cont_reb_lim(cont_reb_lim >=1) = 1;
cont_reb_lim(cont_reb_lim <=0) = 0;
cont_reb_lim(isnan(cont_reb_lim)) = 0;
cont_reb_lim = cont_reb_lim/max(cont_reb_lim(:));
% combine_all(reb_disinh_map==0 & reb_combined==1) = ...
%     combine_all(reb_disinh_map==0 & reb_combined==1) + cont_reb_lim(reb_disinh_map==0 & reb_combined==1);
% combine_all(reb_af_CaKO==0) = combine_all(reb_af_CaKO==0) + cont_reb_lim(reb_af_CaKO==0);
% 
% % Wrong area marked
% % combine_all(combine_all>0 & combine_all<1) = NaN;
% 
% % Spontaneous + rebound
% combine_all(reb_disinh_map==3 & reb_combined==1) = 0;
% INDS = combine_all==3;
% combine_all(INDS) = combine_all(INDS) + cont_reb_lim(INDS);
% 
% combine_all(reb_disinh_map==2 & reb_combined==0) = 4 + reg_af_bef(reb_disinh_map==2 & reb_combined==0)/100;
% combine_all(reb_disinh_map==2 & reb_combined==2) = 4 + reg_af_bef(reb_disinh_map==2 & reb_combined==2)/100;
% combine_all(reb_disinh_map==3 & reb_combined==0) = 2 + spontaneous(reb_disinh_map==3 & reb_combined==0)/100;
% combine_all(reb_disinh_map==3 & reb_combined==2) = 2 + spontaneous(reb_disinh_map==3 & reb_combined==2)/100;

% combine_all(reb_af_CaKO>=ch_th) = 1;

Disinh = 3-spontaneous/100;
Disinh(cont_reb_lim>0)=0;
fr_norm = 1 - firing_rate/max(firing_rate(:));
% combine_all(reb_combined==2)= 0 + fr_norm(reb_combined==2);
combine_all(reb_combined==1)= 0 + cont_reb_lim(reb_combined==1);
% % IND = reb_disinh_map==1 & reb_combined==0;
% % combine_all(IND) = 1 + reb_af_CaKO(IND)/max(reb_af_CaKO(:));
IND = reb_disinh_map==3 & reb_combined==0;
% combine_all(IND) = 0 + fr_norm(IND);
% combine_all = Disinh + cont_reb_lim;
IND = combine_all==0 & firing_rate>=0;
% combine_all(IND) = 0 + fr_norm(IND);
combine_all(reb_disinh_map==0 & reb_combined==0) = NaN;

% Rebound referenced

figure;
colormap(RGB)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),combine_all(sel_GS,sel_GC),[-0.05,1])
cRange = caxis; % save the current color range
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
temp_ax = colorbar();
temp_ax.Box = 'off';
temp_ax.TickDirection = 'out';
temp_ax.FontSize = 14;
temp_ax.Label.String = 'Proportion of trials with rebound mode';
% CMAP = colormap(gca);
% CMAP = [0.4,0.4,0.4;CMAP];
% colormap(CMAP)
% temp_ax.Ticks = [-0.1,temp_ax.Ticks];
ind_ticks_fr = find(temp_ax.Ticks>=0 & temp_ax.Ticks<=1);
min_fr = min(firing_rate(:));
max_fr = max(firing_rate(:));
% ticks_frs = linspace(ceil(max_fr),floor(min_fr),length(ind_ticks_fr));
% for idx = 1:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_frs(idx),'%.0f');
% end
% ind_ticks_fr = find(temp_ax.Ticks>=1);
% ticks_prob = linspace(0,1,length(ind_ticks_fr));
% 
% for idx = 2:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_prob(idx),'%.1f');
% end

firing_rate_int = flipud(round(firing_rate/5)*5);
% firing_rate_int = flipud(firing_rate);
firing_rate(firing_rate>=30)=30;
contour_vals = [1:2:7,10:5:30];
hold on
[C,h] = contour(G_CX(sel_GC),G_SNr(sel_GS),firing_rate, contour_vals, 'ShowText', 'off');
clabel(C, h, 'FontSize', 12, 'Color', 'white')
% clabel(C)
h.LineWidth=1;
h.LineColor='white';
caxis(cRange);  % set the color range to the previous one 
title('Proportion of rebound')
fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-rebound-fr-based-1col-isocontour-jet'])
savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-rebound-fr-based-1col-isocontour-jet.fig'])

% Disinhibition referenced

combine_all_dis = (reg_af_bef - cont_reb_lim*100)/100;
combine_all_dis(combine_all_dis<0) = 0;
combine_all_dis(reb_disinh_map==0 & reb_combined==0) = NaN;

figure;
colormap(RGB)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),combine_all_dis(sel_GS,sel_GC),[-0.05,1])
cRange = caxis; % save the current color range
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
temp_ax = colorbar();
temp_ax.Box = 'off';
temp_ax.TickDirection = 'out';
temp_ax.FontSize = 14;
temp_ax.Label.String = 'Proportion of trials with disinhibition mode';
% CMAP = colormap(gca);
% CMAP = [0.4,0.4,0.4;CMAP];
% colormap(CMAP)
% temp_ax.Ticks = [-0.1,temp_ax.Ticks];
ind_ticks_fr = find(temp_ax.Ticks>=0 & temp_ax.Ticks<=1);
min_fr = min(firing_rate(:));
max_fr = max(firing_rate(:));
% ticks_frs = linspace(ceil(max_fr),floor(min_fr),length(ind_ticks_fr));
% for idx = 1:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_frs(idx),'%.0f');
% end
% ind_ticks_fr = find(temp_ax.Ticks>=1);
% ticks_prob = linspace(0,1,length(ind_ticks_fr));
% 
% for idx = 2:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_prob(idx),'%.1f');
% end

firing_rate_int = flipud(round(firing_rate/5)*5);
firing_rate(firing_rate>=30)=30;

hold on
[~,h] = contour(G_CX(sel_GC),G_SNr(sel_GS),firing_rate,contour_vals, 'ShowText', 'off');
h.LineWidth=1;
h.LineColor='white';
clabel(C, h, 'FontSize', 12, 'Color', 'white')
caxis(cRange);  % set the color range to the previous one 
% title('Proportion of disinhibition')
fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-disinh-fr-based-1col-isocontour-jet'])
savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-disinh-fr-based-1col-isocontour-jet.fig'])

% Not disinhibition
combine_all_dis = (reg_af_bef - cont_reb_lim*100)/100;
combine_all_dis(combine_all_dis<0) = 0;
combine_all_dis(reb_disinh_map==0 & reb_combined==0) = NaN;
IND = combine_all>0 | combine_all_dis>0 | isnan(combine_all);
combine_all_dis(~IND) = -0.03;
RGB = [nan_c;0,0,0;Aut];
figure;
colormap(RGB)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),combine_all_dis(sel_GS,sel_GC),[-0.05,1])
cRange = caxis; % save the current color range
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
temp_ax = colorbar();
temp_ax.Box = 'off';
temp_ax.TickDirection = 'out';
temp_ax.FontSize = 14;
temp_ax.Label.String = 'Proportion of trials with disinhibition mode';
% CMAP = colormap(gca);
% CMAP = [0.4,0.4,0.4;CMAP];
% colormap(CMAP)
% temp_ax.Ticks = [-0.1,temp_ax.Ticks];
ind_ticks_fr = find(temp_ax.Ticks>=0 & temp_ax.Ticks<=1);
min_fr = min(firing_rate(:));
max_fr = max(firing_rate(:));
% ticks_frs = linspace(ceil(max_fr),floor(min_fr),length(ind_ticks_fr));
% for idx = 1:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_frs(idx),'%.0f');
% end
% ind_ticks_fr = find(temp_ax.Ticks>=1);
% ticks_prob = linspace(0,1,length(ind_ticks_fr));
% 
% for idx = 2:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_prob(idx),'%.1f');
% end

firing_rate_int = flipud(round(firing_rate/5)*5);
firing_rate(firing_rate>=30)=30;

hold on
[~,h] = contour(G_CX(sel_GC),G_SNr(sel_GS),firing_rate,contour_vals, 'ShowText', 'off');
h.LineWidth=1;
h.LineColor='white';
clabel(C, h, 'FontSize', 12, 'Color', 'white')
caxis(cRange);  % set the color range to the previous one 
% title('Proportion of disinhibition')
fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-disinh-fr-based-2col-isocontour-disinh-jet'])
savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-disinh-fr-based-2col-isocontour-disinh-jet.fig'])

% Not rebound

combine_all(reb_combined==1)= 0 + cont_reb_lim(reb_combined==1);
% % IND = reb_disinh_map==1 & reb_combined==0;
% % combine_all(IND) = 1 + reb_af_CaKO(IND)/max(reb_af_CaKO(:));
IND = reb_disinh_map==3 & reb_combined==0;
% combine_all(IND) = 0 + fr_norm(IND);
% combine_all = Disinh + cont_reb_lim;
IND = combine_all==0 & firing_rate>=0;
% combine_all(IND) = 0 + fr_norm(IND);
combine_all(reb_disinh_map==0 & reb_combined==0) = NaN;

IND = combine_all>0 | combine_all_dis>0 | isnan(combine_all);
combine_all(~IND) = -0.03;
RGB = [nan_c;0,0,0;Aut];
figure;
colormap(RGB)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),combine_all(sel_GS,sel_GC),[-0.05,1])
cRange = caxis; % save the current color range
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
temp_ax = colorbar();
temp_ax.Box = 'off';
temp_ax.TickDirection = 'out';
temp_ax.FontSize = 14;
temp_ax.Label.String = 'Proportion of trials with rebound mode';
temp_ax.Ticks = [0,.5,1];
temp_ax.TickLabels = {'Disinhibition','Dis = Reb','Rebound'};
% CMAP = colormap(gca);
% CMAP = [0.4,0.4,0.4;CMAP];
% colormap(CMAP)
% temp_ax.Ticks = [-0.1,temp_ax.Ticks];
ind_ticks_fr = find(temp_ax.Ticks>=0 & temp_ax.Ticks<=1);
min_fr = min(firing_rate(:));
max_fr = max(firing_rate(:));
% ticks_frs = linspace(ceil(max_fr),floor(min_fr),length(ind_ticks_fr));
% for idx = 1:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_frs(idx),'%.0f');
% end
% ind_ticks_fr = find(temp_ax.Ticks>=1);
% ticks_prob = linspace(0,1,length(ind_ticks_fr));
% 
% for idx = 2:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_prob(idx),'%.1f');
% end

firing_rate_int = flipud(round(firing_rate/5)*5);
firing_rate(firing_rate>=30)=30;

hold on
[~,h] = contour(G_CX(sel_GC),G_SNr(sel_GS),firing_rate,contour_vals, 'ShowText', 'off');
h.LineWidth=1;
h.LineColor='white';
clabel(C, h, 'FontSize', 12, 'Color', 'white')
caxis(cRange);  % set the color range to the previous one 
% title('Proportion of disinhibition')
fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-disinh-fr-based-2col-isocontour-reb-jet'])
savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-disinh-fr-based-2col-isocontour-reb-jet.fig'])

%% Variability of the first spikes

RGB = [nan_c;Aut];
figure;
colormap(RGB)
imagesc(G_CX(sel_GC),G_SNr(sel_GS),var_firstspk(sel_GS,sel_GC),[-0.05,60])
cRange = caxis; % save the current color range
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
temp_ax = colorbar();
temp_ax.Box = 'off';
temp_ax.TickDirection = 'out';
temp_ax.FontSize = 14;
temp_ax.Label.String = [{'Variability of the first thalamic spike'}, {'after the decrease'}];
% temp_ax.Ticks = [0,.5,1];
% temp_ax.TickLabels = {'Disinhibition','Dis = Reb','Rebound'};
% CMAP = colormap(gca);
% CMAP = [0.4,0.4,0.4;CMAP];
% colormap(CMAP)
% temp_ax.Ticks = [-0.1,temp_ax.Ticks];
% ind_ticks_fr = find(temp_ax.Ticks>=0 & temp_ax.Ticks<=1);
% min_fr = min(firing_rate(:));
% max_fr = max(firing_rate(:));
% ticks_frs = linspace(ceil(max_fr),floor(min_fr),length(ind_ticks_fr));
% for idx = 1:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_frs(idx),'%.0f');
% end
% ind_ticks_fr = find(temp_ax.Ticks>=1);
% ticks_prob = linspace(0,1,length(ind_ticks_fr));
% 
% for idx = 2:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_prob(idx),'%.1f');
% end

% firing_rate_int = flipud(round(firing_rate/5)*5);
% firing_rate(firing_rate>=30)=30;

% hold on
% [~,h] = contour(G_CX(sel_GC),G_SNr(sel_GS),firing_rate,contour_vals, 'ShowText', 'off');
% h.LineWidth=1;
% h.LineColor='white';
% clabel(C, h, 'FontSize', 12, 'Color', 'white')
% caxis(cRange);  % set the color range to the previous one 
% title('Proportion of disinhibition')
fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-std-fr-based-2col-isocontour-reb-jet'])
savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-std-fr-based-2col-isocontour-reb-jet.fig'])

%% Median

RGB = [nan_c;Aut];
figure;
colormap(RGB)
% median_firstspk = median_firstspk - 1000;
imagesc(G_CX(sel_GC),G_SNr(sel_GS),median_firstspk(sel_GS,sel_GC),[-0.05,60])
cRange = caxis; % save the current color range
xlabel('Excitatory input strength, g_{CX\rightarrow TC} (nS/\mum^2)')
ylabel('Inhibitory input strength, g_{SNr\rightarrow TC} (nS/\mum^2)')
GCA = gca;
GCA.Box = 'off';
GCA.TickDir = 'out';
GCA.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.XLabel.FontSize = 14;
GCA.XTick = [min(G_CX),GCA.XTick];
GCA.YTick = [min(G_SNr),GCA.YTick];
temp_ax = colorbar();
temp_ax.Box = 'off';
temp_ax.TickDirection = 'out';
temp_ax.FontSize = 14;
temp_ax.Label.String = [{'Median of the first thalamic spike'}, {'after the decrease'}];
% temp_ax.Ticks = [0,.5,1];
% temp_ax.TickLabels = {'Disinhibition','Dis = Reb','Rebound'};
% CMAP = colormap(gca);
% CMAP = [0.4,0.4,0.4;CMAP];
% colormap(CMAP)
% temp_ax.Ticks = [-0.1,temp_ax.Ticks];
% ind_ticks_fr = find(temp_ax.Ticks>=0 & temp_ax.Ticks<=1);
% min_fr = min(firing_rate(:));
% max_fr = max(firing_rate(:));
% ticks_frs = linspace(ceil(max_fr),floor(min_fr),length(ind_ticks_fr));
% for idx = 1:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_frs(idx),'%.0f');
% end
% ind_ticks_fr = find(temp_ax.Ticks>=1);
% ticks_prob = linspace(0,1,length(ind_ticks_fr));
% 
% for idx = 2:length(ind_ticks_fr)
%     temp_ax.TickLabels{ind_ticks_fr(idx)} = num2str(ticks_prob(idx),'%.1f');
% end

% firing_rate_int = flipud(round(firing_rate/5)*5);
% firing_rate(firing_rate>=30)=30;

% hold on
% [~,h] = contour(G_CX(sel_GC),G_SNr(sel_GS),firing_rate,contour_vals, 'ShowText', 'off');
% h.LineWidth=1;
% h.LineColor='white';
% clabel(C, h, 'FontSize', 12, 'Color', 'white')
% caxis(cRange);  % set the color range to the previous one 
% title('Proportion of disinhibition')
fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-median-fr-based-2col-isocontour-reb-jet'])
savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-median-fr-based-2col-isocontour-reb-jet.fig'])