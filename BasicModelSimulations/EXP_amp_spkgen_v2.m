%% Date 06.09.17

% The function generates Poisson spike trains with an exponential amplitude
% distribution. The parameters needed for the spike train are:
% N: Number of spike trains
% f: Freuqency of each
% corr: Degree of correlation
% T_vec: Time vector

%% Date modified 05.01.18

% Considering the fact that the amplitude distribution is a discrete exp.
% decaying distribution the formula in the book chapter by B. Staude
% requires the 1st and 2nd moments of the discrete bounded distribution NOT
% the continuous distribution. So in this version the discrete formula has
% been taken into account.

function spktimes_exp = EXP_amp_spkgen_v2(expmeanrate,N,f,T_vec)


% clear
% close all
% N = 30;             % Population size
% f = 50;             % Population firing rate
% corr = 0.5;           % Average pair-wise correlation among population
% T_vec = 0:0.01:1000;

Ampdist = exppdf(1:N,expmeanrate);
Ampdist = Ampdist/sum(Ampdist);      % Normalizing to have AUC=1

Amps = 1:N;                       % Amplitude which are equal to number of inputs

numindepspks = round((N*f)/sum(Amps.*Ampdist)); % Computing how many independent Poisson spike we need

rndnumexp = numindepspks*Ampdist;

ampsvals = round(rndnumexp);

spktimes_exp = [];
spktrains = [];
spkid = [];

for ampid = 1:length(ampsvals)
    
    if ampsvals(ampid) ~= 0
        
        [spktimes,~] = MIP_imp_v4_beta(0,1,ampsvals(ampid),T_vec);
        spktimes_exp = [spktimes_exp;repmat(spktimes,ampid,1)];
    end
end

% figure;
% % histogram(ceil(Ampdist),[1:1:30])
% title(['Correlation Coefficient = ',num2str(corr)])

% [spktimes,~] = MIP_imp_v4_beta(0,1,round(carrrate),T_vec);
% 
% spktimes_exp = [];
% 
% for amp_id = length(Ampdist):-1:1
%     
%     if length(spktimes) >= round(Ampdist(amp_id)*carrrate)
%         selected_ids = randperm(length(spktimes),round(Ampdist(amp_id)*carrrate));
%         spktimes_exp = [spktimes_exp;repmat(spktimes(selected_ids),amp_id,1)];
%         spktimes(selected_ids) = [];
%     end
% end
    