%% Date created 15.08.18 by M. Mohagheghi

% This script generates the results for the transmission quality for
% exponential and binomial amplitude distributions

% Exponential

% Running the simulation; 1st input arg. specifies which chunk in
% parameter space should be run and 2nd arg. specifies how many chunks
% should the parameter space be devided into

disp(['Running simulations to compute transmission quality', ...
      ' for exponential amplitude distribution ...'])

res_dir = TCmodel_func_bwfor_modEXPgen(100, 560);
vis_res_lumped_mats(res_dir, 'EXP')