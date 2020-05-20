function Plot_Geom_and_Fields(...
    figPos,dx,coeffs,w,w_bndBox,...
    k,z_poles,z_samp,z_samp_weights,z_ast,N_2)

PlotGeometry(figPos,w,z_poles,z_samp,z_ast)

%----------------------------------------------------------------------
% Create grid
n = floor((w_bndBox(2)-w_bndBox(1))/dx); gx = w_bndBox(1) + dx*(0:n); % grids
n = floor((w_bndBox(4)-w_bndBox(3))/dx); gy = w_bndBox(3) + dx*(0:n);    
[X Y] = meshgrid(gx, gy);
Z = X + 1i*Y;

%----------------------------------------------------------------------
% Evaluate on grid
[u_s, di] = EvalGrid_MPS(w,k,coeffs,Z(:),...
    N_2,z_samp_weights,z_poles,z_ast);

u_s = reshape(u_s, size(Z));
di = reshape(di, size(Z));
    
% u_in = u_in_(Z);
% u_tot = u_in+u_s;
u_tot = u_s;

%----------------------------------------------------------------------
% Plot fields
subplot(3,1,2)
PlotField_MPS(gx,gy,real(u_tot),di,w_bndBox)

subplot(3,1,3)
PlotField_MPS(gx,gy,imag(u_tot),di,w_bndBox)

end
