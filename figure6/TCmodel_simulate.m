%% Modified 16.03.17

% Solved issues: more clarity in the discription, smaller output files for
% recording the excitatory and inhibitory spikes, speed enhancement, NEMO
% compatible

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

function dir_name_res = TCmodel_simulate(job_id,num_jobs)

    disp(['THE JOB ID IS: ',num2str(job_id),' OUT OF ',num2str(num_jobs)])

    curr_dir = pwd;
    sub_dir_name = strsplit(curr_dir,'/');
    sub_dir_name = sub_dir_name{end};

%     ws_dir = '/work/ws/nemo/fr_mm1108-Rebound-0';
    ws_dir = pwd;

    dir_name = [ws_dir,'/',sub_dir_name];

    %addpath(pwd,'~/exp-data')
    %% Ranges of variations for both CX and SNr
    dt = 0.01;
    simtime = 1500;
    T = 0:dt:simtime;

    mov_onset = 1000;

    N_CX = 1;
    N_SNr = 30;
    F_CX = 100;
    F_SNr = 50;

    dir_name = [ws_dir,'/',sub_dir_name,'-',num2str(F_CX,'%.0f'),'-detailed/'];

    n_trials_var = length(N_CX);

    % F_CX = 10*ones(size(T));
    % F_SNr = 80*sigmf(T,[-0.1 mov_onset]);

    corr_CX = 0;
    corr_SNr = 0;
    G_SNr_all = 0.1:0.5:4;  %GS
    G_SNr = G_SNr_all;
    G_CX = 0.1:0.5:4;
    num_trials = 1;   %NT
    num_trials_vec = 1:num_trials;     %NT

    % rebound_spk = zeros(length(corr_SNr),length(N_SNr)*length(F_SNr));
    % all_reb_spk = rebound_spk;

    %% Loading data

    % jit_percents = [1:20 50]/100;
    corr_vals = 0%0:0.05:1;   %JV

    % Variable initialization

    NT_GS_JV_TF_all = combvec(G_SNr,corr_vals,G_CX,num_trials_vec);
    NT_GS_JV_TF = NT_GS_JV_TF_all(:,job_id:num_jobs:end);

    spk_bef_mov = zeros(1,size(NT_GS_JV_TF,2));
    spk_af_mov = zeros(1,size(NT_GS_JV_TF,2));
    corrsnr = zeros(1,size(NT_GS_JV_TF,2));
    CaKO_spk_bef_mov = zeros(1,size(NT_GS_JV_TF,2));
    CaKO_spk_af_mov = zeros(1,size(NT_GS_JV_TF,2));

    T_spk_bef_mov       = zeros(1,size(NT_GS_JV_TF,2));
    T_spk_af_mov        = zeros(1,size(NT_GS_JV_TF,2));
    T_CaKO_spk_bef_mov  = zeros(1,size(NT_GS_JV_TF,2));
    T_CaKO_spk_af_mov   = zeros(1,size(NT_GS_JV_TF,2));
    T_reb_spk_wo_exc    = zeros(1,size(NT_GS_JV_TF,2));
    reb_spk_wo_exc      = zeros(1,size(NT_GS_JV_TF,2));
    sub_th_exc_ind      = zeros(1,size(NT_GS_JV_TF,2));

    SNr_spk     = [];
    CX_spk      = [];
    SNr_spkd    = [];
    CX_spkd     = [];
    TH_spk      = [];
    T_KO_TH_spk = [];
    no_inh_TH_spk = [];
    no_exc_TH_spk = [];
    exc_inh_mo_TH_spk = [];


    % comb_trial_num = NT_GS_JV_TF(1,:);
    comb_G_SNr = NT_GS_JV_TF(1,:);
    comb_jit_val = NT_GS_JV_TF(2,:);
    comb_G_CX = NT_GS_JV_TF(3,:);
    comb_trials = NT_GS_JV_TF(4,:);

%     parcheck_dir = [dir_name,'/checkfiles/'];
% 
%     if exist(parcheck_dir,'dir') ~= 7
%         mkdir(parcheck_dir)
%     end
% 
%     if job_id - 1 == 0
%         parpool('local',str2double(getenv('MOAB_PROCCOUNT')))
%         fl_name = ['job-',num2str(job_id),'.mat'];
%         save([parcheck_dir,fl_name],'job_id')
%     else
%         fl_name = ['job-',num2str(job_id-1),'.mat'];
%         while exist([parcheck_dir,fl_name],'file') ~= 2
%             pause(120)
%             disp('Not generated yet!')
%         end
%         system(['rm -f ',parcheck_dir,fl_name])
%         parpool('local',str2double(getenv('MOAB_PROCCOUNT')))
%         fl_name = ['job-',num2str(job_id),'.mat'];
%         save([parcheck_dir,fl_name],'job_id')
%     end

    parfor S = 1:size(NT_GS_JV_TF,2)    % Loop over experimental trials

        disp([comb_G_SNr(S),comb_jit_val(S),comb_G_CX(S),comb_trials(S)])
        [spk_af_mov(S),spk_bef_mov(S),...
        CaKO_spk_af_mov(S),CaKO_spk_bef_mov(S),...
        T_spk_af_mov(S),T_spk_bef_mov(S),...
        T_CaKO_spk_af_mov(S),T_CaKO_spk_bef_mov(S),T_reb_spk_wo_exc(S),...
        sub_th_exc_ind(S),reb_spk_wo_exc(S),...
        cx_spk_time,snr_spk_time,cx_spk_timed,snr_spk_timed,...
        th_spk,th_spk_Ca_KO,th_spk_no_inh,th_spk_excinh_mo,th_spk_no_exc] = ...
            TC_model_CX_SNr_cond_changed_parfor_opt(N_SNr,...
                            F_SNr,comb_jit_val(S),comb_G_SNr(S),...
                            N_CX,F_CX,corr_CX,comb_G_CX(S),...
                            T,mov_onset);
        CX_spk = [CX_spk,[cx_spk_time;...
                          uint32(comb_trials(S)*ones(1,length(cx_spk_time)));...
                          uint32(S*ones(1,length(cx_spk_time)))]];
        SNr_spk = [SNr_spk,[snr_spk_time;...
                            uint32(comb_trials(S)*ones(1,length(snr_spk_time)));...
                            uint32(S*ones(1,length(snr_spk_time)))]];
        CX_spkd = [CX_spkd,[cx_spk_timed;...
                          comb_trials(S)*ones(1,length(cx_spk_time));...
                          S*ones(1,length(cx_spk_time))]];
        SNr_spkd = [SNr_spkd,[snr_spk_timed;...
                            comb_trials(S)*ones(1,length(snr_spk_time));...
                            S*ones(1,length(snr_spk_time))]];
        TH_spk = [TH_spk,[th_spk';...
                          comb_trials(S)*ones(1,length(th_spk));...
                          comb_G_CX(S)*ones(1,length(th_spk));...
                          comb_G_SNr(S)*ones(1,length(th_spk))]];

        T_KO_TH_spk = [T_KO_TH_spk,[th_spk_Ca_KO';...
                                    comb_trials(S)*ones(1,length(th_spk_Ca_KO));...
                                    comb_G_CX(S)*ones(1,length(th_spk_Ca_KO));...
                                    comb_G_SNr(S)*ones(1,length(th_spk_Ca_KO))]];

        no_inh_TH_spk = [no_inh_TH_spk,[th_spk_no_inh';...
                                        comb_trials(S)*ones(1,length(th_spk_no_inh));...
                                        comb_G_CX(S)*ones(1,length(th_spk_no_inh));...
                                        comb_G_SNr(S)*ones(1,length(th_spk_no_inh))]];

        no_exc_TH_spk = [no_exc_TH_spk,[th_spk_no_exc';...
                                        comb_trials(S)*ones(1,length(th_spk_no_exc));...
                                        comb_G_CX(S)*ones(1,length(th_spk_no_exc));...
                                        comb_G_SNr(S)*ones(1,length(th_spk_no_exc))]];

        exc_inh_mo_TH_spk = [exc_inh_mo_TH_spk,[th_spk_excinh_mo';...
                                                comb_trials(S)*ones(1,length(th_spk_excinh_mo));...
                                                comb_G_CX(S)*ones(1,length(th_spk_excinh_mo));...
                                                comb_G_SNr(S)*ones(1,length(th_spk_excinh_mo))]];
    end

    dir_name_res = [dir_name,'/res-for-colorplot/'];

    if exist(dir_name_res,'dir') ~= 7
        mkdir(dir_name_res)
    end

    save([dir_name_res 'gsnr-4-8-gcx-0-8-Ncx-100-MIP-pois-' date '-nSNr-' num2str(N_SNr) '-' num2str(job_id)],...
        'spk_af_mov','spk_bef_mov','CaKO_spk_af_mov','CaKO_spk_bef_mov',...
        'T_spk_af_mov','T_spk_bef_mov','T_CaKO_spk_af_mov','T_CaKO_spk_bef_mov',...
        'sub_th_exc_ind','T_reb_spk_wo_exc',...
        'NT_GS_JV_TF','reb_spk_wo_exc','CX_spk','SNr_spk','CX_spkd','SNr_spkd',...
        'mov_onset','T_KO_TH_spk','TH_spk','no_inh_TH_spk','no_exc_TH_spk',...
        'exc_inh_mo_TH_spk')

end
