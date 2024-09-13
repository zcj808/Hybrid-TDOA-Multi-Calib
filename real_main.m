%% Compute TDOA-M, TDOA-S and Data used to compute DOA using IROS data
clear
clc
exp_num=15;
load('thres.mat');
for i=1:exp_num
    [my_tdoa] =get_tdoa(i,thres(i),"my");
    tdoa_s=my_tdoa';
    [tdoa_m] = get_tdoa(i,thres(i),"su");
    save(sprintf("data/%d/tdoa_s.mat",i),'tdoa_s')
    save(sprintf("data/%d/tdoa_m.mat",i),'tdoa_m')
end
%% Generate g_list to implement comparison. id 6 in audio name corresponds to id 0 marked in mic. arr.
clear
clc
exp_num=15;
g_list=[];
tdoa_err=[];
doa_err=[];
for exp_id=1:exp_num
    g = sound_gen([6;6;6],exp_id);
    plot_g(g,"gt",1)
    [tdoa_s_qua,tdoa_m_qua,err_ang]=measure_detection(g);
    tdoa_err=[tdoa_err;mean([tdoa_m_qua;tdoa_s_qua])];
    doa_err=[doa_err;mean(err_ang)];
    g_list=[g_list,g];
end
save('Dataset/real.mat','g_list');
save('Outputs/tdoa_err.mat','tdoa_err')
save('Outputs/doa_err.mat','doa_err')
%%
clear
clc
load('Dataset/real.mat')
exp_num=15;
init_my_err=[];
my_err=[];
for i=1:exp_num
    ori_g=g_list(i);
    total_t=tic;
    while toc(total_t)<ori_g.lim_t
        g=init_estimator(ori_g);
        g.label="final";
        [g,norm_dk,value_f] = GN_Solver(g);
        if norm_dk<g.dk_p || value_f<g.f_p
            init_my_err=[init_my_err;[g.init_rec,toc(total_t)]];
            my_err=[my_err;[g.rec,toc(total_t)]];
            break
        end
    end
end
real_zhang_rec=[my_err(:,1),init_my_err(:,3:7),my_err(:,3:7)];
save("Outputs/real_zhang_rec.mat","real_zhang_rec")
%%
clear
clc
load("Outputs/real_wang_rec.mat")
load("Outputs/real_zhang_rec.mat")
load('Outputs/tdoa_err.mat')
load('Outputs/doa_err.mat')

disp("Real world experiment results")
disp("Ours")
Ours=mean(real_zhang_rec(:,2:end))
disp("Wang's")
Wang=mean(real_wang_rec(:,2:end))
disp("Measurement Errors")
tdoa_err=mean(tdoa_err)
doa_err=mean(doa_err)