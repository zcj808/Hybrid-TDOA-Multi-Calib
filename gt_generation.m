% the true values configuration (Assume frame established by first mic. arr.)
function [g] = gt_generation(tdoa_sigma,doa_sigma,traj_id,arr_M)
    %     arr_M=[2,3];
    fs=16000;
    off=0.1; % maximum of time offset and unit is second
    max_n=1e-4;
    M=size(arr_M,2);
    m=8;  % parameters number of each microphone
    mic_lie_a=[];
    for i=1:M
        mic_lie_ai=-1+2*rand(1,3);
        ang=-pi+2*pi*rand; % -pi~pi
        mic_lie_ai=ang*mic_lie_ai/norm(mic_lie_ai);
%         mic_lie_ai=[0,0,0];
        mic_lie_a=[mic_lie_a;mic_lie_ai];
    end
    if traj_id==1
        x_gt=[3*rand(M,3),mic_lie_a,-off+2*off*rand(M,1),-max_n+2*max_n*rand(M,1)]; % traj1
    elseif traj_id==2
        x_gt=[2*rand(M,1),6*rand(M,1),2*rand(M,1),mic_lie_a,-off+2*off*rand(M,1),-max_n+2*max_n*rand(M,1)]; % traj2
    elseif traj_id==3
        x_gt=[-2+4*rand(M,1),4*rand(M,1),2*rand(M,1),mic_lie_a,-off+2*off*rand(M,1),-max_n+2*max_n*rand(M,1)]; % traj3
    end
    x_gt(1,1:6)=0;

    n=3; % parameters number of each sound    
    % traj1
    gt_loc=3*[0,1,1,0,0,0,1,1;
            0,0,1,1,1,0,0,1;
            0,0,0,0,1,1,1,1];
    K=8; % number of sound sources
    if traj_id==2
        % traj2
        gt_loc=[0,2,2,2,2,2,2,2,2,0;
                0,0,2,4,6,6,4,2,0,0;
                0,0,0,0,0,2,2,2,2,2];
        K=10;
    elseif traj_id==3
        gt_loc=2*[0,1,1,1,0,-1,-1,-1,-1,-1,-1,0,1,1;
                0,0,1,2,2,2,1,0,0,1,2,2,2,1;
                0,0,0,0,0,0,0,0,1,1,1,1,1,1];
        K=14;
    end
    g.x_gt = [x_gt;[gt_loc'+0.5,zeros(K,5)]];

    g.x_gt(1,7)=0; % time offset is relative and clock drift is absoulte
    dt=[0;10*rand(K-1,1)]; % my method starts from 2, su starts from 1
    g.x=zeros(size(g.x_gt));

    % generate TDOA for adjacent sound events
    g.tdoa_s = zeros(K-1,M);
    g.tdoa_sigma=tdoa_sigma;
    cc = 340;  % sound speed，unit is m/s
    for i=1:M
        mic_loc = g.x_gt(i,1:3);
        dri = g.x_gt(i,8);
        for j=2:K
            s_loc = g.x_gt(M+j,1:3);
            pre_s_loc = g.x_gt(M+j-1,1:3);
            g.tdoa_s(j-1,i) = (norm(mic_loc-s_loc)-norm(mic_loc-pre_s_loc))/cc+(1+dri)*dt(j)+normrnd(0,g.tdoa_sigma);
        end
    end

    % generate DOA between mic. arr. and sound.
    g.doa_sigma=doa_sigma; % unit: degree
    g.doa = zeros(K,M,3);
    for j=1:K
        s_loc=g.x_gt(M+j,1:3);
        for i=1:M
            mic_loc=g.x_gt(i,1:3);
            lie_a=g.x_gt(i,4:6);
            R=rotvec2mat3d(lie_a);
            dir_vec=(mic_loc'-s_loc')/norm(mic_loc-s_loc);
            
            err_lie_a=-1+2*rand(1,3);
            err_ang=normrnd(0,g.doa_sigma);
            if err_ang>180
                err_ang=180;
            elseif err_ang<-180
                err_ang=-180;
            end
            err_ang=err_ang/180*pi;
            err_lie_a=err_ang*err_lie_a/norm(err_lie_a);
            err_R=rotvec2mat3d(err_lie_a);
            g.doa(j,i,1:3)=(err_R*R'*dir_vec)';
        end
    end
    
    % generate odometry measurement
    g.odo_sigma = 3e-2;
    g.odo=[];
    for j=2:K
        g.odo = [g.odo,g.x_gt(M+j,1:3)'-g.x_gt(M+j-1,1:3)'+normrnd(0,g.odo_sigma,3,1)];
    end
    
    % generate noise convariance matrix(CM)
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
    
    % other configuration
    g.arr_M=arr_M;g.cc = cc;g.fs=fs;g.dt=dt;g.M = M;g.m = m;g.K = K;g.n = n;
    g.dk_p=1e-2;g.f_p=150;

    % generate su tdoa
    sg=g;
    sg.x_gt(2:sg.M,8)=sg.x_gt(2:g.M,8)-sg.x_gt(1,8);
    sg.x_gt(1,8)=0;
    g.g.tdoa_m=zeros(sg.M-1,sg.K);
    mic1_loc=sg.x_gt(1,1:3);
    for j=1:sg.K
        s_loc = sg.x_gt(g.M+j,1:3);
        for i=2:sg.M
            mic_loc = sg.x_gt(i,1:3);
            off=sg.x_gt(i,7); % assume it's relative to mic1
            dri=sg.x_gt(i,8); % assume it's relative to mic1
            % sound and all mics open simultaneously and after delta T emitting
            % a sound and so on...
            g.tdoa_m(i-1,j)=(norm(mic_loc-s_loc)-norm(mic1_loc-s_loc))/sg.cc+off+sum(sg.dt(1:j))*dri+normrnd(0,sg.tdoa_sigma);
        end
    end
end