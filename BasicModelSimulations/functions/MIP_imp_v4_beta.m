function [spike_times_vec,spike_times] = ...
    MIP_imp_v4_beta(des_corr,N,F,t_vec)

%% This function has been written to generate Poissonian spike train with the
%% the specified correlation coefficient using Multi Interaction Process
%% This is the version modified on 20 Nov. 15 which considers the inconsistencies 
%% between the desired and actual correlation coeffiecients


% Modified on 15.02.16

% The purpose of the modification is that the algorithm had a problem when
% the time vector exceeded 1 sec. the problem was due to the usage of sum
% function to determine the frequency of the signal. In this modification
% the sum function is normalized by the length of the time vector to avoid
% the miscalculation of the frequency.

% des_corr = 0.0; % desired correlation
% N = 100 ; % number of inputs
% F = 50; % Frequency
    if des_corr == 0
        moth_fr = F*N*20;
    else
        moth_fr = F/des_corr;
    end
    F = max(F);
%     t_vec = 0:.1:1000;
    mother_spk_train = spkgen(t_vec,1,moth_fr,0);
    
%     if des_corr <= 0.7
%         while sum(mother_spk_train) < moth_fr
%             mother_spk_train = spkgen(t_vec,1,moth_fr,0);
%         end
%     else
    if des_corr ~= 0
        while sum(mother_spk_train)/max(t_vec)*1000 < floor(moth_fr) || ...
                sum(mother_spk_train)/max(t_vec)*1000 > ceil(moth_fr)
            mother_spk_train = spkgen(t_vec,1,moth_fr,0);
        end
    end
%     end
    
%     while des_corr == 1 && sum(mother_spk_train)~=F
%         mother_spk_train = spkgen_v2(t_vec,1,moth_fr,0);
%     end
%     AA = sum(moth_spk_times);
    F = F * max(t_vec)/1000;
    moth_spk_times = t_vec(mother_spk_train > .5); 
    ind_of_corr = zeros(F,N);
%     disp(F)
%     disp(N)
%     disp(length(moth_spk_times))
%     disp(des_corr)
    for ind_int = 1:N
        ind_of_corr(:,ind_int) = randperm(length(moth_spk_times),F);
    end

    spike_times = moth_spk_times(ind_of_corr);
    spike_times = spike_times';
%     spike_trains = zeros(length(t_vec),N);
%     for ind_trains = 1:size(spike_times,2);
%         [spike_trains(2:end,ind_trains),~] = histcounts(spike_times(:,ind_trains),t_vec);
% %     spike_trains(2:end,ind_trains) = hist_temp.Values;
%     end
%    R2 = R2_measure(spike_times);
    spike_times_vec = reshape(spike_times,numel(spike_times),1);
%     corr_coef_mat = corrcoef(spike_trains);
%     mean_corr_val = mean2(corr_coef_mat(logical(tril(ones(size(corr_coef_mat))))...
%                             &~isnan(corr_coef_mat)));
%     spike_times = reshape(spike_times,[numel(spike_times),1]);
