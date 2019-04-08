%% Date created 30.04.18 by M. Mohagheghi

% This script concatenates most of the variables in result files

function [res_file] = data_concatenate(data_path,res_file)
    
%     if exist(res_file,'dir') ~= 7
%         mkdir(res_file)
%     end
%     data_path = [pwd,'/res-for-colorplot/'];
    mat_files = what(data_path);
    mat_files = mat_files.mat;
    NT_GS_JV_TF = [];
    CaKO_spk_af_mov = [];       % All spikes excluding those caused by Ca-Ch after the decrease
    CaKO_spk_bef_mov = [];      % All spikes excluding those caused by Ca-Ch before the decrease
    reb_spk_wo_exc = [];        % Rebound spike caused by pure inhibition
    spk_af_mov = [];            % All spikes, either rebound or regular after the decrease
    spk_bef_mov = [];           % All spikes, either rebound or regular before the decrease
    sub_th_exc_ind = [];        % Subthreshould excitatory inputs (without inhibition)
    T_CaKO_spk_af_mov = [];
    T_CaKO_spk_bef_mov = [];
    T_reb_spk_wo_exc = [];
    T_spk_af_mov = [];
    T_spk_bef_mov = [];
    fl_inds = [];
    TH_spks = [];
    TH_spks_Ca_KO = [];
    TH_spks_no_inh = [];
    TH_spks_no_exc = [];
    TH_spks_mo_exc_inh = [];


%     res_dir = [pwd,'/colorplots-simple-measures/'];
%     if exist(res_dir,'dir') ~= 7
%         mkdir(res_dir)
%     end

    if exist(res_file,'file') ~= 2

        for filename_ind = 1:length(mat_files)
        %     clear NT_GS_JV_TF G_SNr all_reb_spk rebound_spk actual_cc
            disp(['Reading file: ',mat_files{filename_ind}])
            temp_stuc(filename_ind) = load(fullfile(data_path,mat_files{filename_ind}));
        %     G_SNr = [G_SNr,temp_stuc(filename_ind).G_SNr];
            NT_GS_JV_TF = [NT_GS_JV_TF,temp_stuc(filename_ind).NT_GS_JV_TF];
            spk_af_mov = [spk_af_mov,temp_stuc(filename_ind).spk_af_mov];
            spk_bef_mov = [spk_bef_mov,temp_stuc(filename_ind).spk_bef_mov];
            CaKO_spk_af_mov = [CaKO_spk_af_mov,temp_stuc(filename_ind).CaKO_spk_af_mov];
            CaKO_spk_bef_mov = [CaKO_spk_bef_mov,temp_stuc(filename_ind).CaKO_spk_bef_mov];
            T_spk_af_mov = [T_spk_af_mov,temp_stuc(filename_ind).T_spk_af_mov];
            T_spk_bef_mov = [T_spk_bef_mov,temp_stuc(filename_ind).T_spk_bef_mov];
            T_CaKO_spk_af_mov = [T_CaKO_spk_af_mov,temp_stuc(filename_ind).T_CaKO_spk_af_mov];
            T_CaKO_spk_bef_mov = [T_CaKO_spk_bef_mov,temp_stuc(filename_ind).T_CaKO_spk_bef_mov];
            sub_th_exc_ind = [sub_th_exc_ind,temp_stuc(filename_ind).sub_th_exc_ind];
            reb_spk_wo_exc = [reb_spk_wo_exc,temp_stuc(filename_ind).reb_spk_wo_exc];
            T_reb_spk_wo_exc = [T_reb_spk_wo_exc,temp_stuc(filename_ind).T_reb_spk_wo_exc];
            fl_inds = [fl_inds,ones(size(temp_stuc(filename_ind).spk_af_mov))*filename_ind];
            TH_spks = [TH_spks,temp_stuc(filename_ind).TH_spk];
            TH_spks_Ca_KO = [TH_spks_Ca_KO,temp_stuc(filename_ind).T_KO_TH_spk];
            TH_spks_no_inh = [TH_spks_no_inh,temp_stuc(filename_ind).no_inh_TH_spk];
            TH_spks_no_exc = [TH_spks_no_exc,temp_stuc(filename_ind).no_exc_TH_spk];
            TH_spks_mo_exc_inh = [TH_spks_mo_exc_inh,temp_stuc(filename_ind).exc_inh_mo_TH_spk];
            mov_onset = temp_stuc(filename_ind).mov_onset;
        %     corrsnr = [corrsnr;temp_stuc(filename_ind).corrsnr];
            clear temp_stuc
        end

        save(res_file,'NT_GS_JV_TF','spk_af_mov','spk_bef_mov',...
            'CaKO_spk_af_mov','CaKO_spk_bef_mov','T_spk_af_mov','T_spk_bef_mov',...
            'T_CaKO_spk_af_mov','T_CaKO_spk_bef_mov','sub_th_exc_ind','reb_spk_wo_exc',...
            'T_reb_spk_wo_exc','fl_inds','mat_files',...
            'TH_spks','TH_spks_Ca_KO','TH_spks_no_inh',...
            'TH_spks_no_exc','TH_spks_mo_exc_inh','mov_onset')
    else
        disp(['Data ',res_file,' already exists!'])
    end

end
