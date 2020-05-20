function PlotField_MPS_CQ(gx,gy,u,di,w_bndBox,time_N)




figure
colorbar
for t=1:time_N
    clf
    hold on, grid on
    imagesc(gx,gy,u(:,:,t),'alphadata', ~isnan(di))
    shading interp
%     GoodCAxis(u);
    colormap(jet)
    colorbar
    rotate3d on
    xlim([w_bndBox(1),w_bndBox(2)])
    ylim([w_bndBox(3),w_bndBox(4)])
    axis equal
    axis off;
    
    
%     plot(u(:,t),'b')
% %         plot(u(:,t),'b+')
%     hold on, grid on
% 
%     for j=1:length(cornerInds)
%         plot([cornerInds(j),cornerInds(j)],[-maxYLim,maxYLim],'r')
%     end
%     ylim([-maxYLim,maxYLim])
% 
%     timeStr = sprintf('Time: %.2f', t*dt);
%     title(timeStr);
% %         title(sprintf(strcat(titleStr, timeStr())));
    pause(0.01)
end







hold on, grid on
imagesc(gx,gy,u, 'alphadata', ~isnan(di))
shading interp
GoodCAxis(u);
colormap(jet)
colorbar
rotate3d on
xlim([w_bndBox(1),w_bndBox(2)])
ylim([w_bndBox(3),w_bndBox(4)])
axis equal
axis off;
end