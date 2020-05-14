function PlotGeometry(figPos,w,z_poles,z_samp,z_ast)

hFig = figure;
set(hFig, 'Position', figPos)
subplot(3,1,1)
hold on, grid on
scatter(real(w),imag(w),'k*');
scatter(real(z_poles(:)),imag(z_poles(:)),'r.');
scatter(real(z_samp),imag(z_samp),'g.');
scatter(real(z_ast),imag(z_ast),'b*');
axis equal

end

