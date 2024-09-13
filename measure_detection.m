function [tdoa_s_qua,tdoa_m_qua,err_ang] = measure_detection(g)
    % estimate standard deviation of TDOA-M noise
    tdoa_m_det=zeros(g.M-1,g.K);
    A=[];
    for j=1:g.K
        s_loc = g.x_gt(g.M+j,1:3);
        A=[A;[1,sum(g.dt(1:j))]];
        for i=2:g.M
            mic_loc = g.x_gt(i,1:3);
            tdoa_m_det(i-1,j)=g.tdoa_m(i-1,j)-(norm(mic_loc-s_loc)-norm(s_loc))/g.cc;
        end
    end
    tdoa_m_qua=[];
    noise_est=[];
    for i=1:size(tdoa_m_det,1)
        b=tdoa_m_det(i,:)';
        x=(A'*A)\A'*b;
        tau=x(1);
        delta=x(2);
        for j=1:g.K
            noise_est=[noise_est,tdoa_m_det(i,j)-tau-sum(g.dt(1:j))*delta];
        end
        tdoa_m_qua=[tdoa_m_qua;abs(std(noise_est))];
    end

    % estimate standard deviation of TDOA-S noise
    tdoa_s_det=zeros(g.M,g.K-1);
    for j=2:g.K
        for i=1:g.M
            pre_dx = g.x_gt(i,1:3)-g.x_gt(g.M+j-1,1:3);
            dx = g.x_gt(i,1:3)-g.x_gt(g.M+j,1:3);
            tdoa_s_det(i,j-1)=g.tdoa_s(j-1,i)-g.dt(j)-(norm(dx)-norm(pre_dx))/g.cc;
        end
    end
    tdoa_s_qua=[];
    for i=1:size(tdoa_s_det,1)
        dri_est=sum(tdoa_s_det(i,:))/sum(g.dt);
        new_sample=tdoa_s_det(i,:)-dri_est*g.dt(2:end)';
        tdoa_s_qua=[tdoa_s_qua;abs(std(new_sample))];
    end
    
    err_ang=[];
    % doa noise estimation
    for i=1:g.M
        for j=1:g.K
            mea_doa=squeeze(g.doa(j,i,:));
            mic_loc=g.x_gt(i,1:3);
            s_loc=g.x_gt(g.M+j,1:3);
            R_gt=rotvec2mat3d(g.x_gt(i,4:6));
            if i==1
                gt_doa=R_gt'*(-s_loc'/norm(s_loc));
            else
                gt_doa=R_gt'*(mic_loc'-s_loc')/norm(mic_loc-s_loc);
            end
            % gt_doa
            % mea_doa
            acos_err=min(max(gt_doa'*mea_doa,-1),1);
            % acos(acos_err)*180/pi
            err_ang=[err_ang;acos(acos_err)*180/pi];
        end
    end
end