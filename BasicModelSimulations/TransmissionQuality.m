%% Date created 15.08.18 by M. Mohagheghi

% This script generates the results for the transmission quality for
% exponential and binomial amplitude distributions

addpath('./functions')

% Exponential

% Running the simulation; 1st input arg. specifies which chunk in
% parameter space should be run and 2nd arg. specifies how many chunks
% should the parameter space be devided into

disp(['Running simulations to compute transmission quality', ...
      ' for exponential amplitude distribution ...'])

res_dir_exp = TCmodel_func_bwfor_modEXPgen(1, 1);
vis_res_lumped_mats(res_dir_exp, 'EXP')


disp(['Running simulations to compute transmission quality', ...
      ' for binomial amplitude distribution ...'])

res_dir_mip = TCmodel_func_bwfor(100, 420);
vis_res_lumped_mats(res_dir_mip, 'MIP')

disp(['Showing the summary results ...'])
res_summary(res_dir_exp, res_dir_mip)