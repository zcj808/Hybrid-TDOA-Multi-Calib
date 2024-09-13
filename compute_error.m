% Mean Square Error of microphone positions, time offsets, clock drift
% rates and sound event locations
function [mic_err,euler_err,off_err,dri_err,s_err] = compute_error(g)
    euler_err=[];
    for i=2:g.M
        p=[1;1;1];
        p=p/norm(p);
        R_gt=rotvec2mat3d(g.x_gt(i,4:6));
        R_est=rotvec2mat3d(g.x(i,4:6));
        acos_err=(R_gt*p)'*(R_est*p);
        acos_err=min(max(acos_err,-1),1);
        euler_err=[euler_err,180/pi*acos(acos_err)];
    end
    euler_err=sqrt(mean(euler_err.^2));
    mic_err=g.x(2:g.M,1:3)'-g.x_gt(2:g.M,1:3)';
    mic_err=sum(mic_err.^2);
    mic_err=sqrt(mean(mic_err));

    off_err=sqrt(mean((g.x(2:g.M,7)-g.x_gt(2:g.M,7)).^2));
    dri_err=sqrt(mean((g.x(2:g.M,8)-g.x(1,8)-(g.x_gt(2:g.M,8)-g.x_gt(1,8))).^2));
    s_err=g.x(g.M+1:end,1:3)'-g.x_gt(g.M+1:end,1:3)';
    s_err=sum(s_err.^2);
    s_err=sqrt(mean(s_err));
end