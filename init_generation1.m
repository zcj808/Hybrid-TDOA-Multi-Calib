function [g] = init_generation1(g)
    traj_id=g.traj_id;
    % mic. arr. loc. and s. loc. generate in the calibration space randomly
    if traj_id==1
        g.x(:,1:3)=3*rand(g.M+g.K,3);  % traj1
    elseif traj_id==2
        g.x(:,1:3)=[2*rand(g.M+g.K,1),6*rand(g.M+g.K,1),2*rand(g.M+g.K,1)];  % traj2
    elseif traj_id==3
        g.x(:,1:3)=[-2+4*rand(g.M+g.K,1),4*rand(g.M+g.K,1),2*rand(g.M+g.K,1)];  % traj3
    elseif traj_id==0 % real-world traj.
        g.x(:,1:3)=[-1+2*rand(g.M+g.K,1),-2+4*rand(g.M+g.K,1),-0.5*rand(g.M+g.K,1)];
    end
    % mic. arr. ori. generate completely randomly.
    for i=2:g.M
        mic_lie_ai=-1+2*rand(1,3);
        ang=-pi+2*pi*rand; % -pi~pi
        g.x(i,4:6)=ang*mic_lie_ai/norm(mic_lie_ai);
    end
    g.x(1,1:6)=0;
    g.x(1:g.M,7:8)=zeros(g.M,2);
end

