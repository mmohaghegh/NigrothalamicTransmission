%% TC_model_CX_SNr_cond_changed.m description
% This function is written to simulate the neuron with Poissonian inputs
% both from CX and SNr

% the function is written in a way that some input parameters can be changed.
% Also the output of this function is a binary which tells us whether
% rebound happened and count the number of spike. The algorithm is based on
% minimum value of Ica current. The number of spikes are check in a
% predefined window of size 50 ms which should be corrected.

%% added on 01.07.15

% All other event-unrelated rebound spikes before the event is also
% outputed.

%% Function definition

function [spk_af_mov,spk_bef_mov,...
          CaKO_spk_af_mov,CaKO_spk_bef_mov,...
          T_spk_af_mov,T_spk_bef_mov,...
          T_CaKO_spk_af_mov,T_CaKO_spk_bef_mov,T_reb_spk_wo_exc,...
          sub_th_exc_inp,reb_spk_wo_exc,...
          spk_CX,spk_SNr,spk_CXd,spk_SNrd,...
          th_spk_times,T_CaKO_th_spk_times,...
	      no_inh_th_spk_times,th_spks_times_excinh_mo,...
          th_spk_times_wo_exc] = ...
    TC_model_CX_SNr_cond_changed_parfor_opt...
                                          (n_snr,...
                                          f_snr,...
                                          corr_snr,...
                                          g_snr,...
                                          N_CX,...
                                          F_CX,...
                                          corr_CX,...
                                          g_cx,...
                                          T,...
                                          mov_onset)

%% Global variables

%global dir_name

%% input parameters determination

% switch nargin
%     case 2
%         N_CX = n_cx;
%         N_SNr = n_snr;
%         F_CX = 2;
%         F_SNr = 2;
%         corr_CX = 0;
%         corr_SNr = 0;
%     case 4
%         N_CX = n_cx;
%         N_SNr = n_snr;
%         F_CX = f_cx;
%         F_SNr = f_snr;
%         corr_CX = 0;
%         corr_SNr = 0;
%     case 6
%         N_CX = n_cx;
        N_SNr = n_snr;
%         F_CX = f_cx*ones(size(T));
        F_SNr = f_snr;%*sigmf(T,[-0.1 mov_onset]);
%         corr_CX = corr_cx;
        corr_SNr = corr_snr;
% end

% disp([corr_snr,g_snr])

%% Simulation parameters

stepsize = T(2)-T(1);
simtime = T(end);

%% TC parameters


asg = 200;
bsg = 0.4;
itc = 6;
shi = -80;
shi2 = -90;
dur = 5;
dur2 = 10;
period = 25;
gnabar = 3;
gkbar = 5;
glbar = 0.05;
ena = 50;
ek = -90;
eleak = -70;
gtbar = 5;
qht = 5.5;
tadj = 1;
apr = 4;
apt = 0.3;

%% Initial Conditions

vt0 = 0;
ht0 = 0;
htt0 = 0;

%% Synaptic parameters (SNr to TH)

vsynSNrtoTH = -85;
gsynSNrtoTH = g_snr;    % It is multiplied by 8 because in the origiran model
                    % that is used to investigate the effect of DBS, 8 GPi
                    % neurons provide inputs for a single TC cell
alphaSNrtoTH = 1;
betaSNrtoTH = 0.08;
s0SNrtoTH = 0;

%% Synaptic parameters (SNr to TH)

vsynCXtoTH = 0;
gsynCXtoTH = g_cx;    % It is multiplied by 8 because in the origiran model
                    % that is used to investigate the effect of DBS, 8 GPi
                    % neurons provide inputs for a single TC cell
alphaCXtoTH = 1;
betaCXtoTH = 0.19;%0.5;%0.19;
s0CXtoTH = 0;

%% Constant input to the model
% 
% constant_input = 0;
% 
% step_time = 500;
% %simtime = 2000;
% 
% period = simtime/2;
% pulse_width = 20;
% phase_delay = 500;

%% Poissonian input to the TC model


[spk_SNr,spike_times_SNr] = MIP_imp_v4_beta(corr_SNr,N_SNr,F_SNr,...
                                                T(T<=mov_onset));
[spk_CX,spike_times_CX] = MIP_imp_v4_beta(corr_CX,N_CX,F_CX,...
                                                T);

%% Showing the results

% Pure excitation to understand how does the model neuron respond to
% the excitatory inputs

spk_times_SNr = -1;
spk_times_CX = spk_CX;

gnabar = 3;
gtbar = 5;

simopt = simset('SrcWorkspace','current');
sim('TC_RT2004_SNr',[],simopt)

% Checking whether the excitation itself can lead to TC activation

sub_th_exc_inp = isempty(findpeaks(vth.signals.values,'MINPEAKHEIGHT',-40));

if ~sub_th_exc_inp
    [~,spk_inds] = findpeaks(vth.signals.values,'MINPEAKHEIGHT',-40);
    no_inh_th_spk_times = vth.time(spk_inds);
else
    no_inh_th_spk_times = NaN;
end

% Pure inhibition to understand how does the model neuron respond to
% the inhibitory inputs
gnabar = 3;
gtbar = 5;

spk_times_SNr = spk_SNr;
spk_times_CX = -1;

simopt = simset('SrcWorkspace','current');
sim('TC_RT2004_SNr',[],simopt)


% Chechking whether there is rebound spike with pure inhibition

[~,reb_spk_wo_exc] = findpeaks(vth.signals.values,'MINPEAKHEIGHT',-40);
if isempty(reb_spk_wo_exc)
    T_reb_spk_wo_exc = NaN;
    th_spk_times_wo_exc = NaN;
else
    T_reb_spk_wo_exc = vth.time(reb_spk_wo_exc(1));
    th_spk_times_wo_exc = vth.time(reb_spk_wo_exc)
end
reb_spk_wo_exc = length(reb_spk_wo_exc);


% Excitation and inhibition stops at the mov_onset time to understand 
% the deviations in average firing rate from the scenario where excitation
% is ongoing.

gnabar = 3;
gtbar = 5;

spk_times_SNr = spk_SNr;
spk_times_CX = spk_CX(spk_CX<=mov_onset);

simopt = simset('SrcWorkspace','current');
sim('TC_RT2004_SNr',[],simopt)

[~,spk_inds] = findpeaks(vth.signals.values,'MINPEAKHEIGHT',-40);
if isempty(reb_spk_wo_exc)
    th_spks_times_excinh_mo = NaN;
else
    th_spks_times_excinh_mo = vth.time(spk_inds);
end

% Excitatory and inhibitory inputs to understand the role of excitation

spk_times_CX = spk_CX;
spk_times_SNr = spk_SNr;

gnabar = 3;
gtbar = 5;

sim('TC_RT2004_SNr',[],simopt)

spk_bef_mov = length(findpeaks(vth.signals.values(vth.time<mov_onset),...
    'MINPEAKHEIGHT',-40)); % 26 ms is roghly the latency of the rebound spike
TI = vth.time(vth.time<mov_onset);
[~,IND_SPK] = findpeaks(vth.signals.values(vth.time<mov_onset),...
    'MINPEAKHEIGHT',-40);
if isempty(IND_SPK)
    T_spk_bef_mov = NaN;
else
    T_spk_bef_mov = TI(IND_SPK(1));
end

spk_af_mov = length(findpeaks(vth.signals.values(vth.time>=mov_onset),...
    'MINPEAKHEIGHT',-40));

TI = vth.time(vth.time>=mov_onset);
[~,IND_SPK] = findpeaks(vth.signals.values(vth.time>=mov_onset),...
    'MINPEAKHEIGHT',-40);
if isempty(IND_SPK)
    T_spk_af_mov = NaN;
else
    T_spk_af_mov = TI(IND_SPK(1));
end

% Recording the spike times of the thalamocortical neuron.

[~,spk_inds] = findpeaks(vth.signals.values,'MINPEAKHEIGHT',-40);
if ~isempty(spk_inds)
    th_spk_times = vth.time(spk_inds);
else
    th_spk_times = NaN;
end



% T-type Ca Knocked-out maintaining the same excitability

gtbar = 0;
gnabar = 6;

sim('TC_RT2004_SNr',[],simopt)

CaKO_spk_bef_mov = length(findpeaks(vth.signals.values(vth.time<mov_onset),...
    'MINPEAKHEIGHT',-40)); % 26 ms is roghly the latency of the rebound spike

TI = vth.time(vth.time<mov_onset);
[~,IND_SPK] = findpeaks(vth.signals.values(vth.time<mov_onset),...
    'MINPEAKHEIGHT',-40);
if isempty(IND_SPK)
    T_CaKO_spk_bef_mov = NaN;
else
    T_CaKO_spk_bef_mov = TI(IND_SPK(1));
end

CaKO_spk_af_mov = length(findpeaks(vth.signals.values(vth.time>=mov_onset),...
    'MINPEAKHEIGHT',-40));
TI = vth.time(vth.time>=mov_onset);
[~,IND_SPK] = findpeaks(vth.signals.values(vth.time>=mov_onset),...
    'MINPEAKHEIGHT',-40);
if isempty(IND_SPK)
    T_CaKO_spk_af_mov = NaN;
else
    T_CaKO_spk_af_mov = TI(IND_SPK(1));
end

% Recording the spike times of the thalamocortical neuron after T KO.
[~,spk_inds] = findpeaks(vth.signals.values,'MINPEAKHEIGHT',-40);
if ~isempty(spk_inds)
    T_CaKO_th_spk_times = vth.time(spk_inds);
else
    T_CaKO_th_spk_times = NaN;
end

spk_CXd = spk_CX';
spk_SNrd = spk_SNr';
spk_CX = uint32(spk_CX'*100);
spk_SNr = uint32(spk_SNr'*100);
