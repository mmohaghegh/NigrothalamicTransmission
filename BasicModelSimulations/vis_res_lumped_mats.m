%% Script for visualizing the results from the unified for loop
% clear
% 
% close all

function [outputres] = vis_res_lumped_mats(res_path, amp_dist_name)
    show_samples = false;
    res_dir = [res_path,'/res-for-colorplot/'];
    mat_files = what(res_dir);
    mat_files = mat_files.mat;
    G_SNr = [];
    NT_GS_JV_TF = [];
    all_reb_spk = [];
    rebound_spk = [];
    actual_cc = [];
    for filename_ind = 1:length(mat_files)
        temp_stuc(filename_ind) = load([res_dir,mat_files{filename_ind}]);
        NT_GS_JV_TF = [NT_GS_JV_TF,temp_stuc(filename_ind).NT_GS_JV_TF];
        all_reb_spk = [all_reb_spk;temp_stuc(filename_ind).all_reb_spk];
        rebound_spk = [rebound_spk;temp_stuc(filename_ind).rebound_spk];
        clear temp_stuc
    end

    G_SNr = sort(unique(NT_GS_JV_TF(1,:)));

    corr_val = unique(NT_GS_JV_TF(2,:));

    ALL_REB = zeros(length(G_SNr),length(corr_val));
    REB_SPK = zeros(length(G_SNr),length(corr_val));
    OVR = zeros(length(G_SNr),length(corr_val));

    rebound_spk(rebound_spk>=1) = 1;
    for G_ind = 1:length(G_SNr)
        for B = 1:length(corr_val)

                ind = find(NT_GS_JV_TF(1,:)==G_SNr(G_ind) & ...
                           abs(NT_GS_JV_TF(2,:)-corr_val(B))<1e-6);
                ALL_REB(G_ind,B) = mean(all_reb_spk(ind,:));
                REB_SPK(G_ind,B) = mean(rebound_spk(ind,:));
                OVR_temp = (rebound_spk(ind,:))./(rebound_spk(ind,:)+all_reb_spk(ind,:));
                OVR_temp(isnan(OVR_temp))=0;
                OVR(G_ind,B) = mean(OVR_temp);

        end
    end
    outputres = fullfile(res_path,['plotting-params-', amp_dist_name]);
    save(outputres, 'OVR', 'G_SNr', 'corr_val')

    figure;
    colormap('jet')
    imagesc(corr_val,G_SNr,OVR)
    set(gca,'FontSize',14)
    box off
    set(gca,'TickDir','out')
    temp_ax = colorbar();
    temp_ax.Box = 'off';
    temp_ax.TickDirection = 'out';
    temp_ax.FontSize = 12;
    % box off

    xlabel('Correlation Coefficient')
    ylabel('G_S_N_r_t_o_T_C (nS/\mum^2)')
    
    outputfig = fullfile(res_path,['TQ-', amp_dist_name]);
    fig_print(gcf, outputfig)

    if show_samples

        figure;
        colormap('jet')
        imagesc(corr_val,G_SNr,OVR)
        set(gca,'FontSize',14)
        box off
        set(gca,'TickDir','out')
        temp_ax = colorbar();
        temp_ax.Box = 'off';
        temp_ax.TickDirection = 'out';
        temp_ax.FontSize = 12;
        % box off

        xlabel('Correlation Coefficient')
        ylabel('G_S_N_r_t_o_T_C (nS/\mum^2)')
        
        fig_print(gcf,'Colorplot-corr')

        %% Visualizing the membrane potential traces

        figure;

        all_tr_dir = [pwd,'/'];
        num_inps = 30;
        tr_dir = [all_tr_dir,'/voltage-traces-and-inputs/'];

        all_G = 0.3;
        trial_num = randperm(100,1)%15;

        % Specific variables 1st
        G = all_G;
        cc = 0;
        col = [0.498,0,0];

        file_name = [num2str(cc*100),...
            '-',num2str(G*100)];
        load([tr_dir,file_name,'.mat'])
        spk_times = mem_v_traces(trial_num).spike_times;
        mov_onset = mem_v_traces(trial_num).mov_onset;

        vth = TC_model_SNr_inps(G,spk_times);

        % figure;
        subplot(2,3,1)
        plot(vth.time-mov_onset,vth.signals.values,'Color',col,'LineWidth',2)
        xlim([-500,100])
        ylim([-100,0])
        hold on
        axis_temp = axis;
        line([0 0],[axis_temp(end-1),axis_temp(end)],'Color','black',...
            'LineStyle','--','LineWidth',1.5)
        box off
        set(gca,'TickDir','out')
        ylabel([{'Membrane'},{'Voltage (mV)'}])
        set(gca,'XTick',[])
        set(gca,'Position',get(gca,'Position')+[0 -.03 0 0])
        set(gca,'FontSize',14)

        subplot(2,3,1+3)
        raster_spk_times(spk_times-mov_onset,col)

        xlim([-500,100])
        ylim([0,num_inps])
        set(gca,'YTick',[0,num_inps])
        hold on
        axis_temp = axis;
        line([0 0],[axis_temp(end-1),axis_temp(end)],'Color','black',...
            'LineStyle','--','LineWidth',1.5)
        box off
        set(gca,'TickDir','out')

        xlabel('Time (ms)')
        ylabel('Input ID')

        % Specific variables 2nd
        G = 0.3;
        cc = 0.4;
        col = [.18,1,.81];

        file_name = [num2str(cc*100),...
            '-',num2str(G*100)];
        load([tr_dir,file_name,'.mat'])

        for tr_ind = 1:100
            allrebs(tr_ind) = mem_v_traces(tr_ind).all_reb;
        end

        tr_inds_sel = find(allrebs>1);
        trial_num = tr_inds_sel(randperm(length(tr_inds_sel),1));

        spk_times = mem_v_traces(trial_num).spike_times;
        mov_onset = mem_v_traces(trial_num).mov_onset;



        vth = TC_model_SNr_inps(G,spk_times);
        set(gca,'FontSize',14)

        % figure;
        subplot(2,3,2)
        plot(vth.time-mov_onset,vth.signals.values,'Color',col,'LineWidth',2)
        xlim([-500,100])
        ylim([-100,0])
        hold on
        axis_temp = axis;
        line([0 0],[axis_temp(end-1),axis_temp(end)],'Color','black',...
            'LineStyle','--','LineWidth',1.5)
        box off
        set(gca,'TickDir','out')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'Position',get(gca,'Position')+[0 -.03 0 0])
        set(gca,'FontSize',14)
        subplot(2,3,2+3)
        raster_spk_times(spk_times-mov_onset,col)
        xlim([-500,100])
        ylim([0,num_inps])
        set(gca,'YTick',[])
        hold on
        axis_temp = axis;
        line([0 0],[axis_temp(end-1),axis_temp(end)],'Color','black',...
            'LineStyle','--','LineWidth',1.5)
        box off
        set(gca,'TickDir','out')

        xlabel('Time (ms)')

        % Specific variables 3rd
        G = 0.3;
        cc = 0.7;
        col = [0,0,0.87];

        file_name = [num2str(cc*100),...
            '-',num2str(G*100)];
        load([tr_dir,file_name,'.mat'])


        for tr_ind = 1:100
            allrebs(tr_ind) = mem_v_traces(tr_ind).all_reb;
        end

        tr_inds_sel = find(allrebs>5);
        trial_num = tr_inds_sel(randperm(length(tr_inds_sel),1));

        spk_times = mem_v_traces(trial_num).spike_times;
        mov_onset = mem_v_traces(trial_num).mov_onset;

        vth = TC_model_SNr_inps(G,spk_times);
        set(gca,'FontSize',14)

        % figure;
        subplot(2,3,3)
        plot(vth.time-mov_onset,vth.signals.values,'Color',col,'LineWidth',2)
        xlim([-500,100])
        ylim([-100,0])
        hold on
        axis_temp = axis;
        line([0 0],[axis_temp(end-1),axis_temp(end)],'Color','black',...
            'LineStyle','--','LineWidth',1.5)
        box off
        set(gca,'TickDir','out')
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'Position',get(gca,'Position')+[0 -.03 0 0])
        set(gca,'FontSize',14)
        subplot(2,3,3+3)
        raster_spk_times(spk_times-mov_onset,col)
        xlim([-500,100])
        ylim([0,num_inps])
        set(gca,'YTick',[])
        hold on
        axis_temp = axis;
        line([0 0],[axis_temp(end-1),axis_temp(end)],'Color','black',...
            'LineStyle','--','LineWidth',1.5)
        box off
        set(gca,'TickDir','out')

        xlabel('Time (ms)')
        set(gca,'FontSize',14)
        % ylabel('Input ID')
        fig_print(gcf,'memtraces')
    end
end