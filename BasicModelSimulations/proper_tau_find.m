%% Date Created 05.01.18 by M. Mohagheghi

% This function gets a vector consists of average correlation coefficients
% which are desired and then find the proper exponential parameters, tau 
% (the decay rate) leading to the desired correlation for a specific N
% (number of spike trains).

function [comp_corr_vec,comp_tau_vec] = proper_tau_find(N,corr_res)

    tau_vec = 0:0.00001:10;
    int_cnt = 0;
    
    for tau_ind = 1:length(tau_vec)
        pred_cc(tau_ind) = moments_exp(N,tau_vec(tau_ind));
    end
    
    des_corr_vec = 0:corr_res:1;
    
    for d_c_ind = 1:length(des_corr_vec)
%         if des_corr_vec(d_c_ind) > 0
%             temp_ind = find(pred_cc >= (des_corr_vec(d_c_ind)-corr_res/10) & ...
%                             pred_cc <= (des_corr_vec(d_c_ind)+corr_res/10),1);
        [diff_val,temp_ind] = min(abs(pred_cc - des_corr_vec(d_c_ind)));
%         else
%             temp_ind = find(pred_cc <= (des_corr_vec(d_c_ind)+corr_res/10),1);
%         end
        if diff_val <= corr_res/10
            int_cnt = int_cnt + 1;
            comp_corr_vec(int_cnt) = des_corr_vec(d_c_ind);
            comp_tau_vec(int_cnt) = tau_vec(temp_ind);
        end
    end
end