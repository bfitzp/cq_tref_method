function Plot_BdyField_CQ(P,Z,cornerInds,u,T_Min,time_N,dt)

    [~,wi] = intersect(Z,P);

    figure
    colorbar    
    maxYLim = 1.5*max(u(:));
    for t=1:time_N
        clf
        plot(u(:,t),'b')
        hold on, grid on
        
        for j=1:length(wi)
            plot([wi(j),wi(j)],[-maxYLim,maxYLim],'r')
        end
        ylim([-maxYLim,maxYLim])
%         ylim([-0.0574,0.0574])
        
        timeStr = sprintf('Time: %.2f', t*dt);
        title(timeStr);
%         title(sprintf(strcat(titleStr, timeStr())));
        pause(0.01)
    end
end