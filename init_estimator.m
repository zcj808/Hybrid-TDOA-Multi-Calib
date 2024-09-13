function [g] = init_estimator(ori_g)
    ori_g.label="init";
    % 1. Estimate locations and timing parameters using hybrid-TDOA and
    % Odometry  measurement
    total_t=tic;
    while toc(total_t)<ori_g.lim_t
        g = init_generation1(ori_g);
        [g,norm_dk,value_f] = GN_Solver(g);
        if norm_dk<g.dk_p || value_f<g.f_p
            break
        end
    end
    
    % 2. Estimate rotations between i-th arr. and first arr using ICP.
    S1=g.x(g.M+1:end,1:3)'; % Sound loc. in first arr. frame
    cen_s1=mean(S1')';
    qS1=S1-cen_s1*ones(1,g.K);
    for i=2:g.M
        Si=[];
        mic_loc=g.x(i,1:3);
        for j=1:g.K
            s_loc=g.x(g.M+j,1:3);
            dis=norm(mic_loc-s_loc);
            Si=[Si,-dis*squeeze(g.doa(j,i,1:3))];
        end
        cen_si=mean(Si')';
        qSi=Si-cen_si*ones(1,g.K);
        W=zeros(3);
        for j=1:g.K
            W=W+qS1(:,j)*qSi(:,j)';
        end
        [U,S,V] = svd(W);
        R=U*V';
        if det(R)<0
            R=-R;
        end
        g.x(i,4:6)=rotmat2vec3d(R)';
%         pre_err=norm(g.x(i,1:3)-g.x_gt(i,1:3))
%         after_err=norm((cen_s1-R*cen_si)'-g.x_gt(i,1:3))
%         g.x(i,1:3)=(cen_s1-R*cen_si)';
    end
end