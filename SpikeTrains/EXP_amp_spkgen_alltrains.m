%% Date 06.09.17

% The function generates Poisson spike trains with an exponential amplitude
% distribution. The parameters needed for the spike train are:
% N: Number of spike trains
% f: Freuqency of each
% corr: Degree of correlation
% T_vec: Time vector

function [spktimes_exp,spike_train] = EXP_amp_spkgen_alltrains(corr,N,f,T_vec)


% clear
% close all
% N = 30;             % Population size
% f = 50;             % Population firing rate
% corr = 0.5;           % Average pair-wise correlation among population
% T_vec = 0:0.01:1000;

expmeanrate = (corr*(N-1)+1)/2;      % Computing the mu parameter of the exponential distribution
% carrrate = f*N/expmeanrate;

Ampdist = exppdf(0:1:(10*N-1),expmeanrate);% Computing the probability distribution
Ampdist = Ampdist/sum(Ampdist);      % Normalizing to have AUC=1

Amps = 1:1:10*N;                       % Amplitude which are equal to number of inputs

numindepspks = round((N*f)/sum(Amps.*Ampdist)); % Computing how many independent Poisson spike we need
% rndnumexp = exprnd(expmeanrate,1,numindepspks);
rndnumexp = numindepspks*Ampdist;
rndnumexp(N) = sum(rndnumexp(N:end));
rndnumexp(N+1:end) = [];

% amps = histogram(rndnumexp,0:1:N); % With the number of independent spikes, drawing random numbers with the distribution
% ampsvals = amps.Values;

ampsvals = round(rndnumexp);

spktimes_exp = [];
spktrains = [];
spkid = [];

for ampid = 1:length(ampsvals)
    
    if ampsvals(ampid) ~= 0
        
        [spktimes,~] = MIP_imp_v4_beta(0,1,ampsvals(ampid),T_vec);
        spktimes_exp = [spktimes_exp;repmat(spktimes,ampid,1)];
        while length(spktimes) > 0
            if ampid == 1 && floor(length(spktimes)/N) ~= 0
                spkid = [spkid,randperm(N,N)];
                spktrains = [spktrains;spktimes(1:N)];
                spktimes(1:N) = [];
            elseif ampid == 1 && floor(length(spktimes)/N) == 0
                spkid = [spkid,randperm(N,length(spktimes))];
                spktrains = [spktrains;spktimes(1:length(spktimes))];
                spktimes = [];
            else
                spkid = [spkid,randperm(N,ampid)];
                spktrains = [spktrains;repmat(spktimes(1),ampid,1)];
                spktimes(1) = [];
            end
        end
    end
end

for spktr_id = 1:N
    spike_train(spktr_id,:) = logical(sum(T_vec==spktrains(spkid==spktr_id)));
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
    