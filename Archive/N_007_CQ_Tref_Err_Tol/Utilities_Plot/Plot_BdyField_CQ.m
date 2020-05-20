function Plot_BdyField_CQ(P,Z,u,time_N,dt)
  
    %----------------------------------------------------------------------
    % Generate the field for frequencies in upper-right quadrant
    %----------------------------------------------------------------------
    u(:,time_N+2-(1:floor(time_N/2))) = conj(u(:,2:floor(time_N/2)+1));
    
    %----------------------------------------------------------------------
    % Performed scaled inverse FFT
    %----------------------------------------------------------------------
    R = eps^(0.5/(time_N+1));
    u = ifft(u,[],2);
    u = bsxfun(@times,u,R.^(-(0:time_N)));

    [~,wi] = intersect(Z,P);

    %----------------------------------------------------------------------
    % Plot at each timestep
    %----------------------------------------------------------------------
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