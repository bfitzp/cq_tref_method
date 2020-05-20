function PlotField(X,Y,u,x_Max,caxisVals)
hold on, grid on
surf(X,Y,u)
shading interp
colormap(jet)
if ~isempty(caxisVals)
    caxis(caxisVals)
end
colorbar
rotate3d on
xlim([-x_Max,x_Max])
ylim([-x_Max,x_Max])
axis equal
end