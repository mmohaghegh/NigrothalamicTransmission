%% Date Created 05.01.18 by M. Mohagheghi

% This function computes the first and second moments of a discrete
% exponential distribution with N (number of spike trains) and tau (decay
% rate)

function pred_avg_cc = moments_exp(N, tau)

    ks = 1:N;   % All amplitudes
    
    E2_E1 = sum((ks.^2).*exp(-tau*ks))/sum(ks.*exp(-tau*ks));
    
    pred_avg_cc = (E2_E1 - 1)/(N - 1);
    
end