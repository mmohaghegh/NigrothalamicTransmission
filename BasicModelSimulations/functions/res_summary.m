%% Modified on 170915
% Purpose of modification is to additionally plot the threshold of
% correlation at which transmission quality is not good. The threshold now
% is 0.95 meaning that in at least 95% of simulated trials there is only
% specific rebound spike.

%% Script for visualizing the results from the unified for loop
function [] = res_summary(EXP_dir, MIP_dir)

    GSN_sel = 0.70;

    EXP = load(fullfile(EXP_dir, 'plotting-params-EXP.mat'));
    MIP = load(fullfile(MIP_dir, 'plotting-params-MIP.mat'));
    % MIPEXP = load('plotting-params-MIP8EXP2-EXPCC25.mat');
    % EXP = load('plotting-params-EXP.mat');

    figure;
    colormap('jet')
    imagesc(MIP.corr_val,MIP.G_SNr,MIP.OVR)
    GCA = gca;
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.GridColor = [0,0,0];

    box off
    set(gca,'TickDir','out')
    temp_ax = colorbar();
    temp_ax.Box = 'off';
    temp_ax.TickDirection = 'out';
    temp_ax.FontSize = 14;
    temp_ax.Label.String = 'Transmission Quality, TQ';
    % box off

    xlabel('Correlation Coefficient')
    ylabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')


    fig_print(gcf,['Colorplot-corr-MIP-',date])

    figure;
    colormap('jet')
    imagesc(EXP.corr_val,EXP.G_SNr(1:20),EXP.OVR(1:20,:))
    GCA = gca;
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.GridColor = [0,0,0];

    box off
    set(gca,'TickDir','out')
    temp_ax = colorbar();
    temp_ax.Box = 'off';
    temp_ax.TickDirection = 'out';
    temp_ax.FontSize = 14;
    temp_ax.Label.String = 'Transmission Quality, TQ';
    % box off

    xlabel('Correlation Coefficient')
    ylabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')
    fig_print(gcf,['Colorplot-corr-EXP-',date])

%     figure;
%     colormap('jet')
%     imagesc(MIPEXP.corr_val,MIPEXP.G_SNr(1:20),MIPEXP.OVR(1:20,:))
%     GCA = gca;
%     GCA.FontSize = 14;
%     GCA.YLabel.FontSize = 14;
%     GCA.XLabel.FontSize = 14;
%     GCA.GridColor = [0,0,0];
% 
%     box off
%     set(gca,'TickDir','out')
%     temp_ax = colorbar();
%     temp_ax.Box = 'off';
%     temp_ax.TickDirection = 'out';
%     temp_ax.FontSize = 14;
%     temp_ax.Label.String = 'Transmission Quality, TQ';
%     % box off
% 
%     xlabel('Correlation Coefficient')
%     ylabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')
%     fig_print(gcf,['Colorplot-corr-EXP-',date])

    figure;
    plot(MIP.corr_val,MIP.OVR(abs(MIP.G_SNr-GSN_sel)<1e-4,:),'LineWidth',2,'Color','blue')
    hold on
    plot(EXP.corr_val,EXP.OVR(abs(EXP.G_SNr-GSN_sel)<1e-4,:),'LineWidth',2,'Color','red')
%     plot(MIPEXP.corr_val,MIPEXP.OVR(abs(MIPEXP.G_SNr-GSN_sel)<1e-4,:),'LineWidth',2,'Color','green')
    % set(gca,'FontSize',14)
    % title('(1 / (1 + 2))*1')
    GCA = gca;
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.Title.FontSize = 14;
    GCA.TickDir = 'out';
    GCA.GridColor = [0,0,0];
    box off
    set(gca,'TickDir','out')
    xlabel('Correlation Coefficient')
    ylabel('Transmission Quality (TQ)')
    % title('Transmission Quality, TQ','FontWeight','normal')
%     legend('MIP','EXP','MIP&EXP')
    legend('MIP','EXP')
    fig_print(gcf,['G3-TQ-corr-modmeth-',date])

    %% Finding the threshold of correlation at which transmission quality deteriorates
    %% Threshold = 0.95
    G_SNr = MIP.G_SNr;
    for G_ind = 1:length(G_SNr)
        C_thrtemp = find(MIP.OVR(G_ind,:)>=0.95,1,'last');
        if isempty(C_thrtemp)
            C_thrMIP(G_ind) = nan;
        else
            C_thrMIP(G_ind) = MIP.corr_val(C_thrtemp);
        end
        C_thrtemp = find(EXP.OVR(G_ind,:)>=0.95,1,'last');
        if isempty(C_thrtemp)
            C_thrEXP(G_ind) = nan;
        else
            C_thrEXP(G_ind) = EXP.corr_val(C_thrtemp);
        end
%         C_thrtemp = find(MIPEXP.OVR(G_ind,:)>=0.95,1,'last');
%         if isempty(C_thrtemp)
%             C_thrMIPEXP(G_ind) = nan;
%         else
%             C_thrMIPEXP(G_ind) = MIPEXP.corr_val(C_thrtemp);
%         end
    end

    figure;
    plot(G_SNr,C_thrMIP,'LineWidth',2,'Color','blue')
    hold on
    plot(G_SNr,C_thrEXP,'LineWidth',2,'Color','red')
%     plot(G_SNr,C_thrMIPEXP,'LineWidth',2,'Color','green')
    ylim([0,0.65])
    % set(gca,'FontSize',14)
    % title('(1 / (1 + 2))*1')
    GCA = gca;
    GCA.FontSize = 14;
    GCA.YLabel.FontSize = 14;
    GCA.XLabel.FontSize = 14;
    GCA.Title.FontSize = 14;
    GCA.TickDir = 'out';
    GCA.GridColor = [0,0,0];
    box off
    set(gca,'TickDir','out')
    xlabel('Inhibitory input strength, G_{SNr\rightarrow TC} (nS/\mum^2)')
    xlim([0.3,1])
    ylabel('Threshold correlation for clear transmission')
    % title('Transmission Quality (TQ)','FontWeight','normal')
    legend('MIP','EXP')
%     legend('MIP','EXP','MIP&EXP')
    fig_print(gcf,['Threshold-TQ-corr-mod-',date])
end

