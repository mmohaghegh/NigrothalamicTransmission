%% Date Created 30.04.18 By M. Mohagheghi

% This script computes the average firing rate of TH neurons for each
% parameters combination and the visualizes as plots. The raster plot for
% each parameter combination is also included.

function [putative_rebounds, var_first_spk, median_first_spk] = TH_avgfr_several_conds(res_flname,fig_dir,win_width,overlap,with_fig)

    simtime = 1500;
    col = {'blue','red','green','black','magenta'};
    des_corr = 0;       % The map will be for uncorrelated inhibitory inputs

    if exist(fig_dir) ~= 7
        mkdir(fig_dir)
    end

    load(res_flname)
    
    G_SNr_mat = NT_GS_JV_TF(1,:);
    corr_mat = NT_GS_JV_TF(2,:);
    G_CX_mat = NT_GS_JV_TF(3,:);

    G_SNr = sort(unique(NT_GS_JV_TF(1,:)));
    G_CX = sort(unique(NT_GS_JV_TF(3,:)));
    corr_snr = unique(NT_GS_JV_TF(2,:));
    
    putative_rebounds = zeros(length(G_SNr),length(G_CX));

    for Ginh_ind = 1:length(G_SNr)
        
        sprintf('Plotting the results for Gsnr = %.2f',G_SNr(Ginh_ind))
        
        for Gexc_ind = 1:length(G_CX)
            
            sel_inds = abs(TH_spks_Ca_KO(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks_Ca_KO(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times_KO = TH_spks_Ca_KO(1,sel_inds);
            trials_KO    = TH_spks_Ca_KO(2,sel_inds);
            num_trs_KO = length(unique(trials_KO));

            sel_inds = abs(TH_spks_no_inh(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks_no_inh(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times_noinh = TH_spks_no_inh(1,sel_inds);
            trials_noinh    = TH_spks_no_inh(2,sel_inds);
            num_trs_noinh = length(unique(trials_noinh));
            
            sel_inds = abs(TH_spks_no_exc(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks_no_exc(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times_noexc = TH_spks_no_exc(1,sel_inds);
            trials_noexc    = TH_spks_no_exc(2,sel_inds);
            num_trs_noexc = length(unique(trials_noexc));
            
            sel_inds = abs(TH_spks_mo_exc_inh(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks_mo_exc_inh(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times_mo_excinh = TH_spks_mo_exc_inh(1,sel_inds);
            trials_mo_excinh    = TH_spks_mo_exc_inh(2,sel_inds);
            num_trs_mo_excinh = length(unique(trials_mo_excinh));
            
            sel_inds = abs(TH_spks(4,:) - G_SNr(Ginh_ind))<0.001 & ...
                       abs(TH_spks(3,:) - G_CX(Gexc_ind))<0.001;
            
            spk_times = TH_spks(1,sel_inds);
            trials    = TH_spks(2,sel_inds);
            num_trs = length(unique(trials));
            
            [reb_prob] = rebound_detection(spk_times,spk_times_noexc,spk_times_noinh,mov_onset);
            disp(['rebound probability: ',num2str(reb_prob)])
            putative_rebounds(Ginh_ind,Gexc_ind) = reb_prob;
            
            sel_inds = abs(G_CX_mat - G_CX(Gexc_ind))<0.001 & ...
                       abs(G_SNr_mat - G_SNr(Ginh_ind))<0.001;
            T_spk_af_mov_sel = T_spk_af_mov(sel_inds);
            first_spks = T_spk_af_mov_sel(~isnan(T_spk_af_mov_sel));
            first_spks = first_spks - mov_onset;
%             first_spks = spk_times(spk_times>mov_onset);
%             mean_firstspk = median(first_spks);
%             var_firstspk = var(first_spks);
            if isempty(first_spks)
                var_first_spk(Ginh_ind,Gexc_ind) = NaN;
                median_first_spk(Ginh_ind,Gexc_ind) = NaN;
            elseif length(first_spks) == 1
                var_first_spk(Ginh_ind,Gexc_ind) = NaN;
                median_first_spk(Ginh_ind,Gexc_ind) = median(first_spks);
            else
                var_first_spk(Ginh_ind,Gexc_ind) = std(first_spks);
                median_first_spk(Ginh_ind,Gexc_ind) = median(first_spks);
            end
            
            if with_fig
                figure;
    %             subplot(2,1,1)

                [cnt,t_vec] = PSTH_mov_win_fast(spk_times,win_width,overlap,0,simtime,num_trs,1);
                cnt_gexc(Gexc_ind,:) = cnt;
                pl_h = plot(t_vec,cnt,'LineWidth',2,'color',col{1});
                pl_h.Color(4) = 0.5;

                hold on

                [cnt_KO,t_vec_KO] = PSTH_mov_win_fast(spk_times_KO,win_width,overlap,0,simtime,num_trs_KO,1);
                pl_h = plot(t_vec_KO,cnt_KO,'LineWidth',2,'color',col{2});
                pl_h.Color(4) = 0.5;

                [cnt_noinh,t_vec_noinh] = PSTH_mov_win_fast(spk_times_noinh,win_width,overlap,0,simtime,num_trs_noinh,1);
                pl_h = plot(t_vec_noinh,cnt_noinh,'LineWidth',2,'color',col{3});
                pl_h.Color(4) = 0.5;

                [cnt_noexc,t_vec_noexc] = PSTH_mov_win_fast(spk_times_noexc,win_width,overlap,0,simtime,num_trs_noexc,1);
                pl_h = plot(t_vec_noexc,cnt_noexc,'LineWidth',2,'color',col{4});
                pl_h.Color(4) = 0.5;

                [cnt_moei,t_vec_moei] = PSTH_mov_win_fast(spk_times_mo_excinh,win_width,overlap,0,simtime,num_trs_mo_excinh,1);
                pl_h = plot(t_vec_moei,cnt_moei,'LineWidth',2,'color',col{5});
                pl_h.Color(4) = 0.5;

                title(['G_{SNr\rightarrow TC}=',num2str(G_SNr(Ginh_ind),'%.1f'),...
                       ', G_{CX\rightarrow TC}=',num2str(G_CX(Gexc_ind),'%.1f')])
                legend({'inh+exc','inh+exc-Ca-KO','exc','inh','mov-exc-inh'},'Location','northwest')
                xlabel('Time (ms)')
                ylabel('Firing rate (Hz)')
                GCA = gca;
                GCA.FontSize = 14;
                GCA.TickDir = 'out';
                GCA.Box = 'off';

    %             subplot(2,1,2)
    %             raster_trial(spk_times,trials,col{1})
    %             hold on
    % %             raster_trial(spk_times_KO,trials_KO,col{2})
    %             raster_trial(spk_times_noinh,trials_noinh,col{2})
    %             xlim([0,simtime])
    %             xlabel('Time (ms)')
    %             ylabel('Trial #')
    %             GCA = gca;
    %             GCA.FontSize = 14;
    %             GCA.TickDir = 'out';
    %             GCA.Box = 'off';
                fig_print(gcf,fullfile(fig_dir,['gsn',num2str(G_SNr(Ginh_ind)*10,'%.0f'),'-',...
                                                'gcx',num2str(G_CX(Gexc_ind)*10,'%.0f')]))
                close(gcf)
            end
        end
        if with_fig
            figure;
            imagesc(t_vec,G_CX,cnt_gexc)
            ax = colorbar();
            ax.TickDirection = 'out';
            ax.FontSize = 14;
            ax.Label.String = 'Firing rate (Hz)';
            xlabel('Time (ms)')
            ylabel('G_{CX\rightarrow TC}')
            title(['G_{SNr\rightarrow TC} = ',num2str(G_SNr(Ginh_ind),'%.1f')])
            GCA = gca;
            GCA.FontSize = 14;
            GCA.TickDir = 'out';
            GCA.Box     = 'off';
            fig_print(gcf,fullfile(fig_dir,[num2str(G_SNr(Ginh_ind)*10,'%.0f')]))
            close(gcf)
        end
    end
    
end

function [rebound_p] = rebound_detection(inhexc,inh,exc,mov)
    l_time = max(inh);
    interval = abs(mov-l_time);
    fr_inh = sum(inh>=mov & inh<=l_time)/100*1000;
    fr_exc = sum(exc>=mov & exc<=l_time)/100*1000;
    fr_ei  = sum(inhexc>=mov & inhexc<=l_time)/100*1000;
    rebound_p = (fr_ei - fr_exc)/fr_inh;
end