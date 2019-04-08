%% Date Created 30.04.18 by M. Mohagheghi

% This script handles running simulations and analyzing the generated
% simulation files

% Running this function can take up to few weeks when it runs in single
% core mode!
% all_data = TCmodel_simulate(1, 1);

% all_data = './res-for-colorplots/';
% all_data = '/Volumes/Storage/res-for-colorplot-100Hz';
% all_data = './simdata/detailedTHspktimes/res-for-colorplot-100Hz/';
res_dir = '/home/mohaghegh/PhD/Projects/SNr-input-to-the-thalamus/data-rebdisinh-smth-50';
all_data = '/home/mohaghegh/PhD/Projects/SNr-input-to-the-thalamus/data-rebdisinh-smth-50/res-for-colorplot';
conc_data_path = fullfile(res_dir, 'detailedTHspktimes');
fig_path = fullfile(res_dir, 'PSTH-severalconds/100Hz-20ms/');
fig_path_raster = fullfile(res_dir, 'rasterPSTH/100Hz-20ms/');
fig_path_visins = fullfile(res_dir, 'sample-frs-rasters/100Hz');
fig_disinh = fullfile(res_dir, 'Map/100Hz/');
fig_rebdisinh_var = fullfile(res_dir, 'colormap-100Hz-5ms-exc');

% Check whether the paths exist!
check_dir_exist(conc_data_path)
check_dir_exist(fig_path)
check_dir_exist(fig_path_raster)
check_dir_exist(fig_path_visins)
check_dir_exist(fig_disinh)

disp('Concatenating data to a single file ...')
conc_data_file = data_concatenate(all_data, fullfile(conc_data_path, 'concatenated-data-100Hz.mat'));

disp('Plotting the raster plots and average firing rates ...')
% TH_avgfr_raster_visinstances(conc_data_file,fig_path_raster,20,1)
% TH_avgfr_raster(conc_data_file,fig_path_raster,20,1)

disp('Plotting average firing rate ...')
[putative_reb, var_1stspk, med_1stspk] = TH_avgfr_several_conds(conc_data_file,fig_path,20,1,false);

disp('Plotting reb-disinh and baseline firing rate maps ...')
processed_file = rebound_disinhibition_map_func(conc_data_file,putative_reb,var_1stspk,med_1stspk,0.5,fig_disinh);
maps_combination_modcolormap_simple_contour(processed_file, fig_rebdisinh_var)

function [] = check_dir_exist(path)
    if exist(path, 'dir') ~= 7
        mkdir(path)
    end
end