function PlotField_MPS(gx,gy,u,di,w_bndBox)
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