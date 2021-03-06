%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Frequency domain with square geometry
%----------------------------------------------------------------------
clc
clear all
close all

%----------------------------------------------------------------------
% Includes
%----------------------------------------------------------------------
addpath('Utilities')

%----------------------------------------------------------------------
% Corners
%----------------------------------------------------------------------
w_1 = -1-1i;
w_2 = -1+1i;
w_3 = 1+1i;
w_4 = 1-1i;

%----------------------------------------------------------------------
% Wavenumber
%----------------------------------------------------------------------
k = 8;

%----------------------------------------------------------------------
% Poles
%----------------------------------------------------------------------
sigma = 2;
z_polesCorner_N = 12;
z_polesOffset = 0.0;

[z_poles,N_1,w_vec]...
    = GenPoles(...
    z_polesCorner_N,z_polesOffset,sigma,...
    w_1,w_2,w_3,w_4);

%----------------------------------------------------------------------
% Monomials
%----------------------------------------------------------------------
N_2 = 1 * ceil(z_polesCorner_N/2);
z_ast = 1e-14 + 0i*1;

%----------------------------------------------------------------------
% Sample points
%----------------------------------------------------------------------
z_samp_N = 3*z_polesCorner_N;
[z_samp,M] = GenSamplePoints(z_samp_N,...
    z_poles,z_polesCorner_N,w_vec);

%----------------------------------------------------------------------
% System matrix
%----------------------------------------------------------------------
A = GenSysMatrix(k,M,N_1,N_2,z_samp,z_poles,z_ast);
% save('A_','A')

%----------------------------------------------------------------------
% Right hand side
%----------------------------------------------------------------------
opt = {'ps',1};
u_in_ = SelectIncidentField(k,opt);

b = u_in_(z_samp);

%----------------------------------------------------------------------
% Sample system
%----------------------------------------------------------------------
x = A\b;

%----------------------------------------------------------------------
% Parse coefficients
%----------------------------------------------------------------------
a_re = x(1:N_1);
a_im = x(N_1+1:2*N_1);
b_re = x(2*N_1+1:2*N_1+1+N_2);
b_im = x(2*N_1+1+N_2+1:end);

%----------------------------------------------------------------------
% L_inf boundary error
%----------------------------------------------------------------------
err = norm(A*x-b,inf);
fprintf('||Ax - b||_inf: %.16f\n',err);

%----------------------------------------------------------------------
% Evaluate on grid
%----------------------------------------------------------------------
x_Max = 1;
x_Min = -x_Max;

x_step = 0.025;
x = x_Min:x_step:x_Max;
y = x;
[X,Y] = meshgrid(x);
Z = X + 1i*Y;

u_s = EvalGrid(k,Z,...
    a_re,a_im,b_re,b_im,...
    N_1,N_2,z_poles,z_ast);

% u_in = u_in_(Z);
% u_tot = u_in+u_s;
u_tot = u_s;

%----------------------------------------------------------------------
% Plot geometry and fields
%----------------------------------------------------------------------
hFig = figure;
set(hFig, 'Position', [980 100 760 900])
subplot(3,1,1)
hold on, grid on
scatter(real(w_vec),imag(w_vec),'k*');
scatter(real(z_poles(:)),imag(z_poles(:)),'r.');
scatter(real(z_ast),imag(z_ast),'b*');
scatter(real(z_samp),imag(z_samp),'g+');
axis equal

subplot(3,1,2)
caxisVals = [];
GoodCAxis(real(u_tot));
PlotField(X,Y,real(u_tot),x_Max,caxisVals);

subplot(3,1,3)
caxisVals = [];
GoodCAxis(imag(u_tot));
PlotField(X,Y,imag(u_tot),x_Max,caxisVals);
