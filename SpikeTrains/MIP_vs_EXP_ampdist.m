%% Date created 03.10.17

% The purpose of this script is to visualize the amplitude distribution for
% two Poisson processes with two different higher-order correlations.

clear


corr = 0.3;
N = 30;
fr = 50;
T = 0:0.01:1000;
numtr = 100;
spks_mip = [];
spks_exp = [];
sample_spk_exp = [];

for tr_ind = 1:numtr
    [spks,spk_times] = MIP_imp_v4_beta(corr,N,fr,T);
    spks_mip = [spks_mip,spks];
end

sample_spk_mip = spk_times;

for tr_ind = 1:numtr
    [spks,spk_times] = EXP_amp_spkgen_alltrains(corr,N,fr,T);
    spks_exp = [spks_exp,spks];
end

for spktr_ind = 1:size(spk_times,1)
    spkt = T(spk_times(spktr_ind,:));
    sample_spk_exp = [sample_spk_exp,...
                     [spkt;spktr_ind*ones(size(spkt))]];
end

mip_ampdist = ampdist(reshape(spks_mip,1,numel(spks_mip)));
exp_ampdist = ampdist(reshape(spks_exp,1,numel(spks_exp)));

mip_ampdist_cnt = histcounts(mip_ampdist,0:N);
exp_ampdist_cnt = histcounts(exp_ampdist,0:N);

figure;
subplot(2,2,1)
bar(1:N,mip_ampdist_cnt/sum(mip_ampdist_cnt),'blue')
xlim([1,30])
GCA = gca;
GCA.Box = 'off';
GCA.FontSize = 14;
GCA.TickDir = 'out';
GCA.XLabel.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.Title.FontSize = 14;
GCA.Title.FontWeight = 'normal';

xlabel('Amplitudes')
ylabel('Probability')
title('Multiple interaction process')

subplot(2,2,3)
raster_spk_times(sample_spk_mip,'blue')
GCA = gca;
GCA.Box = 'off';
GCA.FontSize = 14;
GCA.TickDir = 'out';
GCA.XLabel.FontSize = 14;
GCA.YLabel.FontSize = 14;
xlim([0,500])
xlabel('Time (ms)')
ylabel('Input ID')

subplot(2,2,2)
bar(1:N,exp_ampdist_cnt/sum(exp_ampdist_cnt),'red')
xlim([1,30])
GCA = gca;
GCA.Box = 'off';
GCA.FontSize = 14;
GCA.TickDir = 'out';
GCA.XLabel.FontSize = 14;
GCA.YLabel.FontSize = 14;
GCA.Title.FontSize = 14;
GCA.Title.FontWeight = 'normal';
xlabel('Amplitudes')
ylabel('Probability')
title('Exponential decaying')

subplot(2,2,4)
for spk_id = 1:size(sample_spk_exp,2)
    line([sample_spk_exp(1,spk_id) sample_spk_exp(1,spk_id)],...
        [sample_spk_exp(2,spk_id)-1 sample_spk_exp(2,spk_id)],'Color','red')
end
GCA = gca;
GCA.Box = 'off';
GCA.FontSize = 14;
GCA.TickDir = 'out';
GCA.XLabel.FontSize = 14;
GCA.YLabel.FontSize = 14;
xlim([0,500])
xlabel('Time (ms)')
ylabel('Input ID')

fig_print(gcf,[pwd,'/MIP-vs-EXP-ampdist-and-raster-C',num2str(corr*10)])