% Plot trajectory of estimated values and true values of g
function []=plot_g(g,label,rot_label) % label='gt' means visualize only true values
    R=[1,0,0;0,-1,0;0,0,-1];    
    if label=='gt'
        if rot_label==1
            g.x_gt(:,1:3)=(R*g.x_gt(:,1:3)')';
        end
        plot3(g.x_gt(1:g.M,1), g.x_gt(1:g.M,2), g.x_gt(1:g.M,3), 'LineStyle','none','Marker','o', 'MarkerSize', 4,'MarkerEdgeColor','r', 'LineWidth',2);
        hold on;
        plot3(g.x_gt(g.M+1:end,1), g.x_gt(g.M+1:end,2), g.x_gt(g.M+1:end,3), 'LineStyle','none','Marker','x', 'MarkerSize', 4.5,'MarkerEdgeColor','b', 'LineWidth',1);
        legend('Mic. gt.','Sound source gt.');
        fac=0.5;
        for i=1:g.M
            if rot_label==1
                gt_loc=g.x_gt(i,1:3)';
                gt_R=fac*R*rotvec2mat3d(g.x_gt(i,4:6));
            end
            quiver3(gt_loc(1), gt_loc(2), gt_loc(3), gt_R(1,1), gt_R(2,1), gt_R(3,1),'r', 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(gt_loc(1), gt_loc(2), gt_loc(3), gt_R(1,2), gt_R(2,2), gt_R(3,2),'g', 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(gt_loc(1), gt_loc(2), gt_loc(3), gt_R(1,3), gt_R(2,3), gt_R(3,3),'b', 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
        end
    else
        plot3(g.x_gt(1:g.M,1), g.x_gt(1:g.M,2), g.x_gt(1:g.M,3), 'LineStyle','none','Marker','o', 'MarkerSize', 4,'MarkerEdgeColor','r', 'LineWidth',2);
        hold on;
        plot3(g.x(1:g.M,1), g.x(1:g.M,2), g.x(1:g.M,3), 'LineStyle','none','Marker','o', 'MarkerSize', 4,'MarkerEdgeColor','g','MarkerFaceColor','g');
        hold on;
        plot3(g.x_gt(g.M+1:end,1), g.x_gt(g.M+1:end,2), g.x_gt(g.M+1:end,3), 'LineStyle','none','Marker','x', 'MarkerSize', 4.5,'MarkerEdgeColor','b', 'LineWidth',1);
        hold on;
        plot3(g.x(g.M+1:end,1), g.x(g.M+1:end,2), g.x(g.M+1:end,3), 'LineStyle','none','Marker','square', 'MarkerSize', 4.5,'MarkerEdgeColor','c');
        hold on;
        
        % first sound loc. is given true value in order to show accumulated
        % error of odometry
        S=[g.x_gt(g.M+1,1:3)'];
        for j=1:g.K-1
            S=[S,S(1:3,j)+g.odo(1:3,j)];
        end
        plot3(S(1,:), S(2,:), S(3,:), 'LineStyle','none','Marker','square', 'MarkerSize', 4.5,'MarkerEdgeColor','y');
        legend('Mic. gt.','Mic. est.','Sound source gt.','Sound source est.','Odometry');

        fac=0.5;
        for i=1:g.M
            est_loc=g.x(i,1:3);
            est_R=fac*rotvec2mat3d(g.x(i,4:6));
            gt_loc=g.x_gt(i,1:3);
            gt_R=fac*rotvec2mat3d(g.x_gt(i,4:6));
            quiver3(est_loc(1), est_loc(2), est_loc(3), est_R(1,1), est_R(2,1), est_R(3,1),'Color', [1, 0.5, 0.5], 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(est_loc(1), est_loc(2), est_loc(3), est_R(1,2), est_R(2,2), est_R(3,2),'Color', [0.5, 1, 0.5], 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(est_loc(1), est_loc(2), est_loc(3), est_R(1,3), est_R(2,3), est_R(3,3),'Color', [0.5, 0.5, 1], 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(gt_loc(1), gt_loc(2), gt_loc(3), gt_R(1,1), gt_R(2,1), gt_R(3,1),'r', 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(gt_loc(1), gt_loc(2), gt_loc(3), gt_R(1,2), gt_R(2,2), gt_R(3,2),'g', 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
            quiver3(gt_loc(1), gt_loc(2), gt_loc(3), gt_R(1,3), gt_R(2,3), gt_R(3,3),'b', 'LineWidth', 2, 'MaxHeadSize', 2,"DisplayName","");
        end
    end
    xlabel('X (m)');ylabel('Y (m)');zlabel('Z (m)');
    grid on;
    view(30,30);
    hold off;
    grid on; axis equal;
end