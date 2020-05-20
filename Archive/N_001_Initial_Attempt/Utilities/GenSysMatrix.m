%----------------------------------------------------------------------
% Generate the system matrix
%----------------------------------------------------------------------
function A = GenSysMatrix(k,M,N_1,N_2,z_samp,z_poles,z_ast)
A = zeros(M,2*N_1+2*N_2+1);

x = 1:M;
y = 1:N_1;
[X,Y] = ndgrid(x,y);

W = z_samp(X)-z_poles(Y);
A(:,1:N_1) = newman_re_(k,W);
A(:,N_1+1:2*N_1) = newman_im_(k,W);

y = 0:N_2;
[X,Y] = ndgrid(x,y);
W_ast = z_samp(X)-z_ast;
A(:,2*N_1+1:2*N_1+1+N_2) = runge_re_(k,W_ast,Y);

y = 1:N_2;
[X,Y] = ndgrid(x,y);
W_ast = z_samp(X)-z_ast;
A(:,2*N_1+1+N_2+1:end) = runge_im_(k,W_ast,Y);
end