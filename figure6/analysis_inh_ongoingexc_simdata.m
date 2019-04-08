%% Date Created 30.04.18 by M. Mohagheghi

% This script handles running simulations and analyzing the generated
% simulation files

% Running this function can take up to few weeks when it runs in single
% core mode!
% all_data = TCmodel_simulate(1, 1);

% all_data = './res-for-colorplots/';
% all_data = '/Volumes/Storage/res-for-colorplot-100Hz';
% all_data = './simdata/detailedTHspktimes/res-for-colorplot-100Hz/';
conc_data_path = './simdata/detailedTHspktimes';
fig_path = './PSTH-severalconds/100Hz-20ms/';
fig_path_raster = './rasterPSTH/100Hz-20ms/';
fig_path_visins = './sample-frs-rasters/100Hz';
fig_disinh = './Map/100Hz/';

% Check whether the paths exist!
check_dir_exist(conc_data_path)
check_dir_exist(fig_path)
check_dir_exist(fig_path_raster)
check_dir_exist(fig_path_visins)
check_dir_exist(fig_disinh)

disp('Concatenating data to a single file ...')
conc_data_file = data_concatenate(all_data, fullfile(conc_data_path, 'concatenated-data-100Hz.mat'));

disp('Plotting the raster plots and average firing rates ...')
TH_avgfr_raster_visinstances(conc_data_file,fig_path_raster,20,1)
% TH_avgfr_raster(conc_data_path,fig_path_raster,20,1)

disp('Plotting average firing rate ...')
[putative_reb, var_1stspk, med_1stspk] = TH_avgfr_several_conds(conc_data_file,fig_path,20,1,false);

disp('Plotting reb-disinh and baseline firing rate maps ...')
rebound_disinhibition_map_func(conc_data_file,putative_reb,var_1stspk,med_1stspk,0.5,fig_disinh)

function [] = check_dir_exist(path)
    if exist(path, 'dir') ~= 7
        mkdir(path)
    end
end