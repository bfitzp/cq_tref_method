function Plot_BdyField_CQ(z_samp,cornerInds,u,T_Min,time_N,dt)
    figure
    colorbar    
    maxYLim = 1.5*max(u(:,time_N));
    for t=1:time_N
        clf
        plot(u(:,t),'b')
%         plot(u(:,t),'b+')
        hold on, grid on
        
        for j=1:length(cornerInds)
            plot([cornerInds(j),cornerInds(j)],[-maxYLim,maxYLim],'r')
        end
%         ylim([-maxYLim,maxYLim])
        ylim([-0.15,0.15])
        
        timeStr = sprintf('Time: %.2f', t*dt);
        title(timeStr);
%         title(sprintf(strcat(titleStr, timeStr())));
        pause(0.01)
    end
end