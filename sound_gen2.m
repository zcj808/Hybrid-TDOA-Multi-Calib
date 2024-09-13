function [g] = sound_gen2(g) % further change configuration of g
    % add orientation and oringin of mic. arr.
    mic_gt=[];
    Rx=g.x_gt(g.M+7,1:3)-g.x_gt(g.M+6,1:3);
    Rx=Rx'/norm(Rx);
    Rz=[0;0;-1];
    Ry=-cross(Rx,Rz);
    R1=[Rx,Ry,Rz]';
    t1=-R1*g.x_gt(1,1:3)';
    g.x_gt(:,1:3)=(R1*g.x_gt(:,1:3)'+t1*ones(1,g.M+g.K))';
    g.S=R1*g.S;
    g.odo=g.S(:,2:end)-g.S(:,1:end-1);
%     for i=1:g.M
%         mici_loc0=g.x_gt(i,1:3);
%         mici_loc3=g.x_gt(g.M+i,1:3);
%         mici_loc1=g.x_gt(2*g.M+i,1:3);
%         dir_x=mici_loc0-mici_loc3;
%         dir_x=dir_x/norm(dir_x);
%         dir31=mici_loc1-mici_loc3;
%         dir_y=dir31-dir_x*(dir31*dir_x');
%         dir_y=dir_y/norm(dir_y);
%         dir_z=cross(dir_x,dir_y);
%         if i==1
%             mic_gt=[mic_gt;zeros(1,g.m)];
%             R1=[dir_x',dir_y',dir_z']';
%             t1=-R1*mici_loc0';
%         else
%             Ri=[dir_x',dir_y',dir_z'];
%             ti=mici_loc0';
%             mic_gt=[mic_gt;[(R1*ti+t1)',rotmat2vec3d(R1*Ri),0,0]];
%         end
%         sin_value=cross([1;0;0],dir_x');
%         if norm(sin_value)>0
%             sin_value=norm(sin_value);
%         else
%             sin_value=-norm(sin_value);
%         end
%         a=asin(sin_value); % -pi~pi
%         Ra=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
%         mic_lie_ai=rotmat2vec3d(Ra);
%         if i==1
%             mic_gt=[mic_gt;zeros(1,g.m)];
%             R1=Ra';
%             t1=-R1*mici_loc0';
%         else
%             mic_gt=[mic_gt;[(R1*mici_loc0'+t1)',rotmat2vec3d(R1*Ra),0,0]];
%         end
%     end
%     s_locs=(R1*g.x_gt(floor(3*g.M)+1:end,1:3)'+t1*ones(1,g.K))';
%     g.x_gt=[mic_gt;[s_locs,zeros(g.K,5)]];
    g.x=zeros(size(g.x_gt));

%     plot_g(g,"gt",1)
%     for i=1:g.M
%         for j=1:g.K
%             mic_loc=g.x_gt(i,1:3);
%             s_loc=g.x_gt(g.M+j,1:3);
%             R=rotvec2mat3d(g.x_gt(i,4:6));
%             squeeze(g.doa(j,i,:))
%             R'*(mic_loc'-s_loc')/norm(mic_loc-s_loc)
%             disp("-----")
%         end
%     end
end