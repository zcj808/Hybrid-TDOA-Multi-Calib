function [g] = sound_gen(mic_ids,id)
    gt_mic_locs=[];
    % generate real experiment data
    arr_loc_ids=[ones(3,1)*[8,7,3];ones(3,1)*[8,7,6];
        ones(3,1)*[2,4,6];ones(3,1)*[2,5,6];
        ones(3,1)*[2,6,1]]; % [loc. id selected by 1-st mic., loc. id selected by 2-st mic., ...]
    arr_loc_id=arr_loc_ids(id,:);
    % compute microphones locations on abs. frame
    arr_M=[6,6,6];
    M=sum(mic_ids>0,'all');
    m=8;
    p_list=[];
    p4_list=[];
    ex_locs=[-56,0,-48,0;-56,-56,-48,-56;
        0,-55,7,-55;58,-55,65,-55;
        58,0,65,0;-3.5,0,3.5,0;
        25,-27,32,-27;-27,-31,-34,-31]/100;
    hig_list=[2,6,12]/100;
    for i=1:length(arr_loc_id)
        p_list=[p_list,[ex_locs(arr_loc_id(i),1:2)';hig_list(i)]];
        p4_list=[p4_list,[ex_locs(arr_loc_id(i),3:4)';hig_list(i)]];
    end
    a=0.035;
    c=cos(pi/3);
    s=sin(pi/3);
    % array shape
    ori_arr_locs=[0,0,0;
        0,-a*c,-a*s;
        0,-a*(1+c),-a*s;
        0,-2*a,0;
        0,-a*(1+c),a*s;
        0,-a*c,a*s];
    R_list=[];
    for i=1:length(arr_M)
        p1=p_list(:,i);
        p4=p4_list(:,i);
        Ry=p1-p4;
        Ry=Ry/norm(Ry);
        Rx=[0;0;1];
        Rz=cross(Rx,Ry);
        R=[Rx,Ry,Rz];
        R_list=[R_list,R];
    end
    for j=1:length(mic_ids(i,:))
        for i=1:length(arr_M)
            R=R_list(:,3*(i-1)+1:3*i);
            p=p_list(:,i)*ones(1,arr_M(i));
            t_arr_locs=R*ori_arr_locs'+p;
            if mic_ids(i,j)>0
                gt_mic_locs=[gt_mic_locs,t_arr_locs(:,mic_ids(i,j))];
            end
        end
    end
    displace=[];
    K=14; % number of sound sources
    n=3;
    % compute and transform sound locations from abs. frame to s. frame
    gt_s_locs=[-92.5+5.5,-92.5+5.5,-36.5-1,40.5-1,92-5.5,43+1,-37+1;
        0+1,-52+1,-86+5.5,-86+5.5,-26-1,38-5.5,38-5.5];
    gt_s_locs=[[gt_s_locs;80*ones(1,7)],[gt_s_locs;43*ones(1,7)]]/100;
    load(sprintf('displacement/displacement_%d.mat',id));
    S=[0;0;0];
    for i=1:K-1
        S=[S,S(:,i)+displace(i,:)'];
    end
    S(3,1:7)=S(3,1:7)+0.8;
    S(3,8:end)=S(3,8:end)+0.43;
    
    % scatter3(gt_s_locs(1,:),gt_s_locs(2,:),gt_s_locs(3,:));
    % hold on
    % scatter3(gt_mic_locs(1,:),gt_mic_locs(2,:),gt_mic_locs(3,:));
    % legend('gt.s','gt.mic')
    
    [R_mea,t_mea]=any2s(S(:,1),S(:,2),S(:,3));
    g.S=R_mea*S+t_mea*ones(1,K);
    
    [R_gt,t_gt]=any2s(gt_s_locs(:,1),gt_s_locs(:,2),gt_s_locs(:,3));
    g.x_gt=[(R_gt*gt_mic_locs+t_gt*ones(1,M))';(R_gt*gt_s_locs+t_gt*ones(1,K))'];
    g.x_gt=[g.x_gt,zeros(M+K,m-3)];
    
    load('dt.mat')
    g.dt=dt(id,:);
    g.dt=g.dt';

    load(sprintf('data/%d/tdoa_s.mat',id))
    g.tdoa_s=tdoa_s;

    load(sprintf('data/%d/tdoa_m.mat',id))
    g.tdoa_m=tdoa_m;

    load(sprintf('data/%d/doa.mat',id))
    g.doa=-double(doa);
    g.doa(:,:,3)=-g.doa(:,:,3);
    g.doa=permute(g.doa, [2, 1, 3]);
    
    % other parameters configuration
    g.M = size(arr_M,2);
    g.m = m;
    g.K = K;
    g.n = n;
    g.tdoa_sigma=1e-4;
    g.doa_sigma=5;
    g.odo_sigma=1e-2;
    % noise convariance matrix
    doa_w=(2-2*cos(g.doa_sigma*pi/180))^0.5;
    tdoaN=M*(K-1)+(M-1)*K;
    doaN=3*M*K;
    odoN=n*(K-1);
    g.W=g.tdoa_sigma^2*eye(tdoaN+doaN+odoN);
    g.W(tdoaN+1:tdoaN+doaN,tdoaN+1:tdoaN+doaN)=doa_w^2*eye(doaN);
    g.W(tdoaN+doaN+1:end,tdoaN+doaN+1:end)=g.odo_sigma^2*eye(odoN);

    % generate noise CM for initial estimation
    W_init=g.tdoa_sigma^2*eye(M*(K-1)+(M-1)*K+3*K+n*(K-1));
    W_init(tdoaN+1:tdoaN+3*K,tdoaN+1:tdoaN+3*K)=doa_w^2*eye(3*K);
    W_init(tdoaN+3*K+1:end,tdoaN+3*K+1:end)=g.odo_sigma^2*eye(n*(K-1));
    g.W_init=W_init;
    g.cc = 340;
    g.dk_p=1e-2;
    g.f_p=150;
    g.traj_id=0;

    % change g twice
    Rx=g.x_gt(g.M+7,1:3)-g.x_gt(g.M+6,1:3);
    Rx=Rx'/norm(Rx);
    Rz=[0;0;-1];
    Ry=-cross(Rx,Rz);
    R1=[Rx,Ry,Rz]';
    t1=-R1*g.x_gt(1,1:3)';
    g.x_gt(:,1:3)=(R1*g.x_gt(:,1:3)'+t1*ones(1,g.M+g.K))';
    g.S=R1*g.S;
    g.odo=g.S(:,2:end)-g.S(:,1:end-1);
    g.x=zeros(size(g.x_gt));
    g.lim_t=2;
    g.id=id;
end