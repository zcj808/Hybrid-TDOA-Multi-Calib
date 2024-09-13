% Compute Jacobian of non-linear least square of TDOA and odometry
function [J,r] = compute_J(g) % Compute Jacobian
    J = [];
    r = [];
    % Jaco. for TDOA-S
    for i=1:g.M
        J_mic = zeros(g.K-1,g.m*g.M-7);
        J_s = zeros(g.K-1,g.n*g.K);
        for j=1:g.K-1
            pre_dx = g.x(i,1:3)-g.x(g.M+j,1:3);
            dx = g.x(i,1:3)-g.x(g.M+j+1,1:3);
            dri = g.x(i,8);
            if i==1
                J_mic(j,1)=[g.dt(j+1)];
            else
                J_mic(j,1+g.m*(i-2)+1:1+g.m*(i-1))=[(dx/norm(dx)-pre_dx/norm(pre_dx))/g.cc,0,0,0,0,g.dt(j+1)];
            end
            J_s(j,g.n*(j-1)+1:g.n*j)=pre_dx/norm(pre_dx)/g.cc;
            J_s(j,g.n*j+1:g.n*j+3)=-dx/norm(dx)/g.cc;
            r = [r;(norm(dx)-norm(pre_dx))/g.cc+(1+dri)*g.dt(j+1)-g.tdoa_s(j,i)];
        end
        J = [J;[J_mic,J_s]];
    end

    % Jaco. for TDOA-M
    for j=1:g.K
        J_mic = zeros(g.M-1,g.m*g.M-7);
        J_s = zeros(g.M-1,g.n*g.K);
        s_loc = g.x(g.M+j,1:3);
        mic_loc1 = g.x(1,1:3);
        dx1=mic_loc1-s_loc;
        dri1 = g.x(1,8);
        J_mic(:,1)=-sum(g.dt(1:j))*ones(g.M-1,1);
        for i=2:g.M
            mic_loc = g.x(i,1:3);
            dx=mic_loc-s_loc;
            off = g.x(i,7);
            dri = g.x(i,8);
            J_mic(i-1,1+g.m*(i-2)+1:1+g.m*(i-1))=[dx/norm(dx)/g.cc,0,0,0,1,sum(g.dt(1:j))]; 
            J_s(i-1,g.n*(j-1)+1:g.n*j)=(dx1/norm(dx1)-dx/norm(dx))/g.cc;
            r = [r;(norm(mic_loc-s_loc)-norm(mic_loc1-s_loc))/g.cc+off+sum(g.dt(1:j))*(dri-dri1)-g.tdoa_m(i-1,j)];
        end
        J = [J;[J_mic,J_s]];
    end

    % Jaco. for DOA
    J_doa=zeros(3*g.M*g.K,g.m*g.M-7+g.n*g.K);
    r_doa=[];
    for j=1:g.K
        s_loc=g.x(g.M+j,1:3);
        for i=1:g.M
            mic_loc=g.x(i,1:3);
            mic_lie_a=g.x(i,4:6);
            R=rotvec2mat3d(mic_lie_a);
            ms_norm=norm(mic_loc-s_loc);
            dir_vec=(mic_loc'-s_loc')/ms_norm;
            Jloc=zeros(3);
            for p=1:3
                for q=p:3
                    if p==q
                        Jloc(p,q)=(1/ms_norm-(mic_loc(p)-s_loc(p))^2/ms_norm^3)/2;
                    else
                        Jloc(p,q)=-(mic_loc(p)-s_loc(p))*(mic_loc(q)-s_loc(q))/ms_norm^3;
                    end
                end
            end
            Jloc=Jloc'+Jloc;
            Jloc=R'*Jloc;
            Jlie=skew(R'*dir_vec);
            if i>1
                J_doa(3*g.M*(j-1)+3*(i-1)+1:3*g.M*(j-1)+3*i,1+g.m*(i-2)+1:1+g.m*(i-2)+3)=Jloc;
                J_doa(3*g.M*(j-1)+3*(i-1)+1:3*g.M*(j-1)+3*i,1+g.m*(i-2)+4:1+g.m*(i-2)+6)=Jlie;
            end
            J_doa(3*g.M*(j-1)+3*(i-1)+1:3*g.M*(j-1)+3*i,g.m*g.M-7+g.n*(j-1)+1:g.m*g.M-7+g.n*j)=-Jloc;
            r_doa=[r_doa;R'*dir_vec-squeeze(g.doa(j,i,1:3))];
        end
    end

    if g.label=="init"
        for j=1:g.K % select doa conerning first arr.
            J=[J;J_doa(3*g.M*(j-1)+1:3*g.M*(j-1)+3,1:end)];
            r=[r;r_doa(3*g.M*(j-1)+1:3*g.M*(j-1)+3)];
        end
    elseif g.label=="final"
        J=[J;J_doa];
        r=[r;r_doa];
    end

    % Jaco. for Odo
    J_odo=zeros(g.n*(g.K-1),g.n*g.K);
    for j=1:g.K-1
        J_odo(g.n*(j-1)+1:g.n*(j-1)+3,g.n*(j-1)+1:g.n*(j+1))=[-eye(3),eye(3)];
        s_now = g.x(g.M+j,1:3)';
        s_next = g.x(g.M+j+1,1:3)';
        r = [r;s_next-s_now-g.odo(:,j)];
    end
    J_odo=[zeros(g.n*(g.K-1),g.m*g.M-7),J_odo];
    J=[J;J_odo];

    if g.label=="init"
        J_init=[J(:,1)];
        for i=1:g.M-1
            J_init=[J_init,J(:,1+g.m*(i-1)+1:1+g.m*(i-1)+3),J(:,1+g.m*(i-1)+7:1+g.m*(i-1)+8)];
        end
        J_init=[J_init,J(:,1+g.m*(g.M-1)+1:end)];
        J=J_init;
    end
end