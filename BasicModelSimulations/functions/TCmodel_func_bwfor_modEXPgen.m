%% Modification date 07.03.16

% Here I modified the script so that it can compute the Pinsky-Rinzel 
% measure of synchrony from the spike trains generated using MIP


%% Date 03.12.15

%% TCmodel_clus_func is function which thoughout the comparison-of-all-3 dir
%% Prepares the intialization like loading experimental data (if necessary)
%% , specifying the parameters required for being sweeped or scanned, doing
%% some analysis on one the experimental data e.g. to extract the trails having
%% almost same rate and same distribution. It also initialize important para-
%% meters for running on clusters.

% The contents of the function mentioned here is to use spike tarins
% based on the MIP mode of spike trains. . To have a
% better results and enough inputs the total simulation time is 2000 ms and
% the movement event occurs at 1500th ms.

% The inputs to this function is job_id and num_jobs which split the
% G_SNr parameter to num_jobs smaller one, so it can be distributed along
% clusters. 'test_bashfile.m' and 'bashfile_gen_bwfor.m' are the scripts be
% relevant for this purpose.
%% 08.07.2015 (For Loop optimization)
% job_id = 1;
% num_jobs = 1;

function [dir_name] = TCmodel_func_bwfor_modEXPgen(job_id,num_jobs)


% The goal of this m-file is to plot the results as Robert wants for his
% reports. In these simulations, there is no cortical inputs and the goal
% is to understand how correlation and conductance of SNr inputs affect
% rebound activity in a thalamocortical cell. Each simulation will be
% repeated 100 times and then as a result, the percentage of reobund
% activity is determined

% addpath(pwd,'~/exp-data')

%% Directory creation

    curr_dir = pwd;
    sub_dir_name = strsplit(curr_dir,'/');
    sub_dir_name = sub_dir_name{end};

    % ws_dir = '/work/ws/nemo/fr_mm1108-Rebound-0';
    ws_dir = pwd;
    % dir_name = [ws_dir,'/',sub_dir_name];
    dir_name = fullfile(ws_dir, 'Results-EXP');

    %% Ranges of variations for both CX and SNr
    dt = 0.01;
    simtime = 1500;
    T = 0:dt:simtime;

    mov_onset = 1000;

    N_CX = 200;
    N_SNr = 30;
    F_CX = 1:0.5:10;
    F_SNr = 50;

    n_trials_var = length(N_CX);

    % F_CX = 10*ones(size(T));
    % F_SNr = 80*sigmf(T,[-0.1 mov_onset]);

    corr_CX = 0;
    corr_SNr = 0;
    G_SNr_all = 0.05:0.05:1;  %GS
    G_SNr = G_SNr_all;
    num_trials = 100;   %NT
    num_trials_vec = 1:num_trials;     %NT

    % rebound_spk = zeros(length(corr_SNr),length(N_SNr)*length(F_SNr));
    % all_reb_spk = rebound_spk;

    %% Loading data

    % jit_percents = [1:20 50]/100;
    % corr_vals = 0:0.05:1;   %JV
    corr_res = 0.05;
    [corr_vals,tau_vec] = proper_tau_find(N_SNr,corr_res);


    NT_GS_JV_TF_all = combvec(G_SNr,corr_vals);
    temp_comb = combvec(G_SNr,tau_vec);
    NT_GS_JV_TF = NT_GS_JV_TF_all(:,job_id:num_jobs:end);

    rebound_spk = zeros(size(NT_GS_JV_TF,2),num_trials);
    all_reb_spk = zeros(size(NT_GS_JV_TF,2),num_trials);
    R2 = zeros(size(NT_GS_JV_TF,2),num_trials);

    comb_G_SNr = NT_GS_JV_TF(1,:);
    comb_jit_val = NT_GS_JV_TF(2,:);
    comb_tau_vec = temp_comb(2,:);

    dir_name_trace = [dir_name,'/voltage-traces-and-inputs/'];
    if exist(dir_name_trace,'dir') ~= 7
        mkdir(dir_name_trace)
    end

    % To run on NEMO cluster
    % parpool('local',str2double(getenv('MOAB_PROCCOUNT')))

    for S = 1:size(NT_GS_JV_TF,2)    % Loop over experimental trials

       [rebound_spk(S,:),all_reb_spk(S,:)] = ...
            TC_model_CX_SNr_cond_changed_parfor_opt_expdistmod(N_SNr,...
                            F_SNr,0,comb_G_SNr(S),...
                            T,mov_onset,comb_jit_val(S),comb_tau_vec(S),...
                            num_trials,dir_name_trace);
    end

    dir_name_cp = [dir_name,'/res-for-colorplot/'];

    if exist(dir_name_cp,'dir') ~= 7
        mkdir(dir_name_cp)
    end

    save([dir_name_cp 'MIP-pois-' date '-nSNr-' num2str(N_SNr) '-' num2str(job_id)],...
        'rebound_spk','all_reb_spk','G_SNr',...
        'num_trials','NT_GS_JV_TF')
end