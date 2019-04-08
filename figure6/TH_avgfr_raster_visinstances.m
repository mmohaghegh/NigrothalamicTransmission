%% Date Created 30.04.18 By M. Mohagheghi

% This script computes the average firing rate of TH neurons for each
% parameters combination and the visualizes as plots. The raster plot for
% each parameter combination is also included.

function [] = TH_avgfr_raster_visinstances(res_flname,fig_dir,win_width,overlap)

    simtime = 1500;
    mov_onset = 1000;
    col = {'blue','red','green'};
    col_cako = [0, 0, 0.5];
    des_corr = 0;       % The map will be for uncorrelated inhibitory inputs

    if exist(fig_dir) ~= 7
        mkdir(fig_dir)
    end

    load(res_flname)
    
    G_SNr_mat = NT_GS_JV_TF(1,:);
    corr_mat = NT_GS_JV_TF(2,:);
    G_CX_mat = NT_GS_JV_TF(3,:);
    
    G_SNr = [2.0, 2.0, 3.4, 1.0];
    G_CX = [0.1, 2.5, 4.0, 1.9];

%     G_SNr = sort(unique(NT_GS_JV_TF(1,:)));
%     G_CX = sort(unique(NT_GS_JV_TF(3,:)));
    corr_snr = unique(NT_GS_JV_TF(2,:));

    for Ginh_ind = 1:length(G_SNr)
        
        sprintf('Plotting the results for Gsnr = %.2f',G_SNr(Ginh_ind))
        
        for Gexc_ind = 1:length(G_CX)
            
            sel_inds = abs(TH_spks(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times = TH_spks(1,sel_inds);
            trials    = TH_spks(2,sel_inds);
            num_trs = length(unique(trials));
            
            sel_inds = abs(TH_spks_Ca_KO(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks_Ca_KO(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times_KO = TH_spks_Ca_KO(1,sel_inds);
            trials_KO    = TH_spks_Ca_KO(2,sel_inds);
            num_trs_KO = length(unique(trials_KO));

%             sel_inds = abs(TH_spks_no_inh(4,:) - G_SNr(Ginh_ind))<0.001 & ...
%                        abs(TH_spks_no_inh(3,:) - G_CX(Gexc_ind))<0.001;
%             
%             spk_times_noinh = TH_spks_no_inh(1,sel_inds);
%             trials_noinh    = TH_spks_no_inh(2,sel_inds);
%             num_trs_noinh = length(unique(trials_noinh));
            
            figure;
            subplot(2,1,1)
            
            [cnt,t_vec] = PSTH_mov_win(spk_times,win_width,overlap,0,simtime,num_trs,1);
            cnt_gexc(Gexc_ind,:) = cnt;
            plot(t_vec-mov_onset,cnt,'LineWidth',2,'color',col{1})
            hold on
            [cnt_KO,t_vec_KO] = PSTH_mov_win(spk_times_KO,win_width,overlap,0,simtime,num_trs_KO,1);
            plot(t_vec_KO-mov_onset,cnt_KO,'LineWidth',2,'color',col_cako)
%             [cnt_noinh,t_vec_noinh] = PSTH_mov_win(spk_times_noinh,win_width,overlap,0,simtime,num_trs_noinh,1);
%             plot(t_vec_noinh,cnt_noinh,'LineWidth',2,'color',col{2})
            title(['G_{SNr\rightarrow TC}=',num2str(G_SNr(Ginh_ind),'%.1f'),...
                   ', G_{CX\rightarrow TC}=',num2str(G_CX(Gexc_ind),'%.1f')])
            legend({'inh+exc','inh+exc-Ca-KO'},'Location','northwest')
            xlabel('Time (ms)')
            ylabel('Firing rate (Hz)')
            ylim([0, 70])
            xlim([-200, 200])
            GCA = gca;
            GCA.FontSize = 14;
            GCA.TickDir = 'out';
            GCA.Box = 'off';
            
            subplot(2,1,2)
            raster_trial(spk_times,trials,col{1})
            hold on
            raster_trial(spk_times_KO,trials_KO,col{2})
%             raster_trial(spk_times_noinh,trials_noinh,col{2})
            xlim([-200,200])
            ylim([0,100])
            xlabel('Time (ms)')
            ylabel('Trial #')
            GCA = gca;
            GCA.FontSize = 14;
            GCA.TickDir = 'out';
            GCA.Box = 'off';
            fig_print(gcf,fullfile(fig_dir,['gsn',num2str(G_SNr(Ginh_ind)*10,'%.0f'),'-',...
                                            'gcx',num2str(G_CX(Gexc_ind)*10,'%.0f')]))
            close(gcf)
        end
%         figure;
%         imagesc(t_vec,G_CX,cnt_gexc)
%         ax = colorbar();
%         ax.TickDirection = 'out';
%         ax.FontSize = 14;
%         ax.Label.String = 'Firing rate (Hz)';
%         xlabel('Time (ms)')
%         ylabel('G_{CX\rightarrow TC}')
%         title(['G_{SNr\rightarrow TC} = ',num2str(G_SNr(Ginh_ind),'%.1f')])
%         GCA = gca;
%         GCA.FontSize = 14;
%         GCA.TickDir = 'out';
%         GCA.Box     = 'off';
%         fig_print(gcf,fullfile(fig_dir,[num2str(G_SNr(Ginh_ind)*10,'%.0f')]))
%         close(gcf)
    end
    
end