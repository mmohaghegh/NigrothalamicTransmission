% Date created: 23.03.17 by Mohammadreza Mohagheghi Nejad

%% This script is going to tell us about how the state-space for rebound 
%% and disinhibition scenario is. The data it uses is the simulation where
%% we apply inhibition having movement-related modulation together with
%% random excitation (f=5Hz, N=10, corr=0)
function res_file = rebound_disinhibition_map_func(res_flname,contaminated_reb,var_firstspk,median_firstspk,reb_th,res_dir)
    load(res_flname)

    sub_th_exc_ind_bin = logical(sub_th_exc_ind);
    spk_bef_mov_bin = logical(spk_bef_mov);
    reb_spk_wo_exc_bin = logical(reb_spk_wo_exc);
    CaKO_spk_af_mov_bin = logical(CaKO_spk_af_mov);
    spk_af_mov_bin = logical(spk_af_mov);

    G_SNr_mat = NT_GS_JV_TF(1,:);
    corr_mat = NT_GS_JV_TF(2,:);
    G_CX_mat = NT_GS_JV_TF(3,:);
    TRs = NT_GS_JV_TF(4,:);

    G_SNr = sort(unique(NT_GS_JV_TF(1,:)));
%     G_SNr = G_SNr(G_SNr<=2);
    G_CX = sort(unique(NT_GS_JV_TF(3,:)));
%     G_CX = G_CX(G_CX<=2);
    corr_snr = unique(NT_GS_JV_TF(2,:));

    ch_th = reb_th*100;

    % res_dir = [pwd,'/colorplots-reb-disinh-map-exc',num2str(max(G_CX)),...
    %             '-inh',num2str(max(G_SNr)),'-Jun17/'];

%     res_dir = [pwd,'/paperfig-180425-200Hz-95percent/'];
    if exist(res_dir,'dir') ~= 7
        mkdir(res_dir)
    end

    des_corr = 0;       % The map will be for uncorrelated inhibitory inputs

    for Ginh_ind = 1:length(G_SNr)
        for Gexc_ind = 1:length(G_CX)

            sel_inds = abs(corr_mat - des_corr)<0.01 & ...
                       abs(G_SNr_mat - G_SNr(Ginh_ind))<0.001 & ...
                       abs(G_CX_mat - G_CX(Gexc_ind))<0.001;

            reg_af_sub(Ginh_ind,Gexc_ind) = sum(~sub_th_exc_ind_bin(sel_inds) & ...
                                            spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds));

            reg_af_reb(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds) & ...
                                            ~reb_spk_wo_exc_bin(sel_inds));

            reg_af_CaKO(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds) & ...
                                            CaKO_spk_af_mov_bin(sel_inds));

            reg_af_bef(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds));

            reb_af_sub(Ginh_ind,Gexc_ind) = sum(sub_th_exc_ind_bin(sel_inds) & ...
                                            spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds));

            reb_af_reb(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds) & ...
                                            reb_spk_wo_exc_bin(sel_inds));

            reb_af_CaKO(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds) & ...
                                            ~CaKO_spk_af_mov_bin(sel_inds));

            reb_af_bef(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                            ~spk_bef_mov_bin(sel_inds));

            no_spk_af_mov(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds));
            
            spontaneous(Ginh_ind,Gexc_ind) = sum(spk_af_mov_bin(sel_inds) & ...
                                                 spk_bef_mov_bin(sel_inds));
            
            th_spk_sel_ind = abs(TH_spks(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                             abs(TH_spks(3,:) - G_CX(Gexc_ind))<0.001;
            
            firing_rate(Ginh_ind,Gexc_ind) = sum(TH_spks(1,th_spk_sel_ind)<=mov_onset)/max(TRs);

        end
    end

    % Contribution of excitation and therefore disinhibition

%     figure;
% 
%     subplot(221)
%     imagesc(G_CX,G_SNr,reg_af_sub)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering subthreshold excitation')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(222)
%     imagesc(G_CX,G_SNr,reg_af_reb)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering rebound inhibition')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(223)
%     imagesc(G_CX,G_SNr,reg_af_CaKO)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering knocked-out Ca channel')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(224)
%     imagesc(G_CX,G_SNr,reg_af_bef)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering only spike before and after')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     set(gcf,'Position',get(0,'ScreenSize'))
%     fig_print(gcf,[res_dir,'corr-',num2str(des_corr*10)])
% 
%     % Contribution of rebound
% 
%     figure;
% 
%     subplot(221)
%     imagesc(G_CX,G_SNr,reb_af_sub)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering subthreshold excitation')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(222)
%     imagesc(G_CX,G_SNr,reb_af_reb)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering rebound inhibition')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(223)
%     imagesc(G_CX,G_SNr,reb_af_CaKO)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering knocked-out Ca channel')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(224)
%     imagesc(G_CX,G_SNr,reb_af_bef)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering only spike before and after')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     set(gcf,'Position',get(0,'ScreenSize'))
%     fig_print(gcf,[res_dir,'reb-corr-',num2str(des_corr*10)])
% 
%     % Contribution of disinhibition, where there is always disinhibitory
%     % mechanism
% 
%     figure;
% 
%     subplot(221)
%     imagesc(G_CX,G_SNr,reg_af_sub>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering subthreshold excitation')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(222)
%     imagesc(G_CX,G_SNr,reg_af_reb>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering rebound inhibition')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(223)
%     imagesc(G_CX,G_SNr,reg_af_CaKO>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering knocked-out Ca channel')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(224)
%     imagesc(G_CX,G_SNr,reg_af_bef>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering only spike before and after')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     set(gcf,'Position',get(0,'ScreenSize'))
%     fig_print(gcf,[res_dir,'95percent-corr-',num2str(des_corr*10)])
% 
%     % Contribution of rebound, where there is always rebound mechanism
% 
%     figure;
% 
%     subplot(221)
%     imagesc(G_CX,G_SNr,reb_af_sub>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering subthreshold excitation')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(222)
%     imagesc(G_CX,G_SNr,reb_af_reb>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering rebound inhibition')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(223)
%     imagesc(G_CX,G_SNr,reb_af_CaKO>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering knocked-out Ca channel')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     subplot(224)
%     imagesc(G_CX,G_SNr,reb_af_bef>=ch_th)
%     xlabel('Excitatory conductance')
%     ylabel('Inhibitory conductance')
%     title('Considering only spike before and after')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
% 
%     set(gcf,'Position',get(0,'ScreenSize'))
%     fig_print(gcf,[res_dir,'95percent-reb-corr-',num2str(des_corr*10)])

    reb_disinh_map = 3*ones(length(G_SNr),length(G_CX));
    reb_disinh_map(reb_af_CaKO>=ch_th)=1;
    reb_disinh_map(logical((reg_af_bef>=ch_th)-(reb_af_CaKO>=ch_th)))=2;
    reb_disinh_map(no_spk_af_mov<ch_th)=0;

    figure;
    colormap('jet')
    imagesc(G_CX,G_SNr,reb_disinh_map)
    xlabel('Excitatory input strength, G_{CX\rightarrow TC} (nS/\mum^2)')
    ylabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')
    GCA = gca;
    GCA.Box = 'off';
    GCA.TickDir = 'out';
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.XTick = [0.05,GCA.XTick];
    GCA.YTick = [0.05,GCA.YTick];
    fig_print(gcf,[res_dir,'rebound-disinhibition-map'])
    
    figure;
    colormap('jet')
    imagesc(G_CX,G_SNr,firing_rate)
    xlabel('Excitatory input strength, G_{CX\rightarrow TC} (nS/\mum^2)')
    ylabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')
    GCA = gca;
    GCA.Box = 'off';
    GCA.TickDir = 'out';
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.XTick = [0.05,GCA.XTick];
    GCA.YTick = [0.05,GCA.YTick];
    temp_ax = colorbar();
    temp_ax.Box = 'off';
    temp_ax.TickDirection = 'out';
    temp_ax.FontSize = 14;
    temp_ax.Label.String = 'Baseline firing rate (Hz)';
    fig_print(gcf,[res_dir,'baseline-FR-map'])
    
    reb_combined = ones(length(G_SNr),length(G_CX));
    reb_combined(isnan(contaminated_reb))=0;
    reb_combined(contaminated_reb>=reb_th)=1;
    reb_combined(contaminated_reb<reb_th)=2;
    
    figure;
    colormap('jet')
    imagesc(G_CX,G_SNr,reb_combined)
    xlabel('Excitatory input strength, G_{CX\rightarrow TC} (nS/\mum^2)')
    ylabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')
    GCA = gca;
    GCA.Box = 'off';
    GCA.TickDir = 'out';
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.XTick = [0.05,GCA.XTick];
    GCA.YTick = [0.05,GCA.YTick];
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
%     temp_ax.FontSize = 14;
%     temp_ax.Label.String = 'Baseline firing rate (Hz)';
    fig_print(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-rebound-fr-based'])
    savefig(gcf,[res_dir,num2str(reb_th*10,'%.0f'),'-rebound-fr-based.fig'])
    
    res_file = fullfile(res_dir,'variables-med-std');
    save(res_file,'reb_combined','reb_disinh_map',...
                                       'contaminated_reb','reb_af_CaKO',...
                                       'reg_af_bef','no_spk_af_mov','spontaneous',...
                                       'contaminated_reb','reb_th',...
                                       'G_CX','G_SNr','firing_rate',...
                                       'var_firstspk','median_firstspk')

%     figure;
%     subplot(121)
%     colormap('jet')
%     imagesc(G_CX,G_SNr,reb_af_CaKO>=ch_th)
%     xlabel('Excitatory conductance (nS/\mum^2)')
%     ylabel('Inhibitory conductance (nS/\mum^2)')
%     title('Rebound plays role')
%     GCA = gca;
%     GCA.Box = 'off';
%     GCA.TickDir = 'out';
%     GCA.FontSize = 12;
% 
%     subplot(122)
%     colormap('jet')
%     imagesc(G_CX,G_SNr,(reg_af_bef>=ch_th)-(reb_af_CaKO>=ch_th))
%     xlabel('Excitatory conductance (nS/\mum^2)')
%     ylabel('Inhibitory conductance (nS/\mum^2)')
%     title('Disinhibition plays role')
%     GCA = gca;
%     GCA.Box = 'off';
%     GCA.TickDir = 'out';
%     GCA.FontSize = 12;
% 
% 
%     set(gcf,'Position',get(0,'ScreenSize'))
%     fig_print(gcf,[res_dir,'Rebound-disinhibition-',num2str(des_corr*10)])
% 
%     % State-Space plot for the rebound-disinhibition map
% 
%     reb_map = reb_af_CaKO>=ch_th;
%     disinh_map = (reg_af_bef>=ch_th)-(reb_af_CaKO>=ch_th);
% 
%     rebmap_x_vals = [];
%     rebmap_y_vals = [];
% 
%     dismap_x_vals = [];
%     dismap_y_vals = [];
% 
%     for row_ind = 1:size(reb_map,1)
% 
%         [~,rcol_ind] = find(reb_map(row_ind,:) == 1,1,'last');
%         [~,dcol_ind] = find(disinh_map(row_ind,:) == 1,1,'first');
%         if ~isempty(rcol_ind)
%             rebmap_x_vals = [rebmap_x_vals,rcol_ind];
%             rebmap_y_vals = [rebmap_y_vals,row_ind];
%         end
%         if ~isempty(dcol_ind)
%             dismap_x_vals = [dismap_x_vals,dcol_ind];
%             dismap_y_vals = [dismap_y_vals,row_ind];
%         end
% 
%     end
%     figure;
%     hold on
%     for ind = 1:numel(rebmap_x_vals)-1
%         line(rebmap_x_vals(ind:ind+1),rebmap_y_vals(ind:ind+1),...
%             'Color','black','LineWidth',2)
%     end
end
