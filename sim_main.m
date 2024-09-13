%% Simulations Set Up
clear
clc7
tdoa_sigmas=[5e-5,1e-4,5e-4]; % TDOA noise sigma
doa_sigmas=[5,10,15];  % unit: degree, 5/10/15
arr_M=[6,6,6,6,6]; % Mic. Number
nums=5;
for doa_i=1:3
    for tdoa_i=1:3
        g_list=[];
        for traj_id=1:3
            for i=1:nums
                ori_g=gt_generation(tdoa_sigmas(tdoa_i),...
                    doa_sigmas(doa_i),traj_id,arr_M);
                ori_g.traj_id=traj_id;
                ori_g.id=nums*(traj_id-1)+i;
                ori_g.lim_t=2; % Time upper limit
                g_list=[g_list,ori_g];
            end
        end
        save(sprintf('Dataset/sim%d.mat',3*(doa_i-1)+tdoa_i),'g_list');
    end
end
%% Zhang's Simulation
clear
clc
for doa_i=1:3
    for tdoa_i=1:3
        index=3*(doa_i-1)+tdoa_i;
        load(sprintf('Dataset/sim%d.mat',index));
        my_err=[];
        init_my_err=[];
        for i=1:size(g_list,2)
            i
            total_t=tic;
            ori_g=g_list(i);
            while toc(total_t)<ori_g.lim_t
                init_g=init_estimator(ori_g);
                g=init_g;
                g.label="final";
                [g,norm_dk,value_f] = GN_Solver(g);
                if norm_dk<g.dk_p || value_f<g.f_p
                    init_my_err=[init_my_err;[g.init_rec,toc(total_t)]];
                    my_err=[my_err;[g.rec,toc(total_t)]];
%                     subplot(2,2,1);
%                     plot_g(init_g,'no');
%                     title('our method init.','FontSize',20)
%                     subplot(2,2,3);
%                     plot_g(g,'no');
%                     title('our method result','FontSize',20)
                    break
                end
            end
        end
        zhang_rec=[my_err(:,1),init_my_err(:,3:7),my_err(:,3:7)];
        save(sprintf("Outputs/zhang_rec%d.mat",index),"zhang_rec")
    end
end
%% Quick Analysis: err=loc. err.(cm), angle err.(deg), off. err.(0.1ms), dri. err.(us)
clear
clc
id=1; % id=3*(doa_id-1)+tdoa_id and doa_id/tdoa_id corresponding to three noise levels of DOA/TDOA
load(sprintf('Outputs/zhang_rec%d.mat',id));
load(sprintf('Outputs/wang_rec%d.mat',id));
z=zhang_rec;
w=wang_rec;
init_zhang_err=[1e2*iqr_mean(z(:,2)),iqr_mean(z(:,3)),1e4*iqr_mean(z(:,4)),...
    1e6*iqr_mean(z(:,5))]
init_wang_err=[1e2*iqr_mean(w(:,2)),iqr_mean(w(:,3)),1e4*iqr_mean(w(:,4)),...
    1e6*iqr_mean(w(:,5))]
zhang_err=[1e2*iqr_mean(z(:,7)),iqr_mean(z(:,8)),1e4*iqr_mean(z(:,9)),...
    1e6*iqr_mean(z(:,10))]
wang_err=[1e2*iqr_mean(w(:,7)),iqr_mean(w(:,8)),1e4*iqr_mean(w(:,9)),...
    1e6*iqr_mean(w(:,10))]

function [iqr_means]=iqr_mean(data)
        q1 = prctile(data, 25);
        q3 = prctile(data, 75);
        iqr_value = q3 - q1;
        lower_bound = q1;
        upper_bound = q3;
        in_iqr = data(data >= lower_bound & data <= upper_bound);
        if ~isempty(in_iqr)
            iqr_means = mean(in_iqr);
        else
            iqr_means = NaN; 
        end
end