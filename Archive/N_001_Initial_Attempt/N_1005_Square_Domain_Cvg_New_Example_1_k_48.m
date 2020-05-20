%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Convergence of error for frequency domain with square geometry
% as number of degrees of freedom increase (poles and monomials)
%----------------------------------------------------------------------
clc
clear all
% close all

%----------------------------------------------------------------------
% Includes
%----------------------------------------------------------------------
addpath('Utilities')

%----------------------------------------------------------------------
% Corners
%----------------------------------------------------------------------
w_1 = -1-1i;
w_2 = -1.0+1.0i;
w_3 = 1+1i;
w_4 = 1-1i;

%----------------------------------------------------------------------
% Wavenumber
%----------------------------------------------------------------------
k = 48 + 0i;
% k = 0 + 25i;

%----------------------------------------------------------------------
% Poles
%----------------------------------------------------------------------
z_polesCorner_Min = 0;
z_polesCorner_Max = 92;
z_polesCorner_Step = 4;
z_polesCorner_Range = ...
    z_polesCorner_Min:z_polesCorner_Step:z_polesCorner_Max;
z_polesCorner_Range_N = length(z_polesCorner_Range);

%----------------------------------------------------------------------
% Iterate with an increasing number of DoFs
%----------------------------------------------------------------------
A_cond = zeros(z_polesCorner_Range_N,1);
A__cond = zeros(z_polesCorner_Range_N,1);
x_norm = zeros(z_polesCorner_Range_N,1);
x_reg_norm = zeros(z_polesCorner_Range_N,1);
err = zeros(z_polesCorner_Range_N,1);
N_vec = zeros(z_polesCorner_Range_N,1);

    
%----------------------------------------------------------------------
% Sample points
w_vec = [w_1;w_2;w_3;w_4];
% z_samp_N = 2*1;
[z_samp,z_samp_weights,M] = GenSamplePoints_New(w_vec);  
    
for j=1:z_polesCorner_Range_N
    %----------------------------------------------------------------------
    % Poles
    sigma = 0.4;
    z_polesCorner_N = z_polesCorner_Range(j);
    z_polesOffset = -0.02;
    
    [z_poles,N_1,w_vec]...
        = GenPoles(...
        z_polesCorner_N,z_polesOffset,sigma,...
        w_1,w_2,w_3,w_4);
    
    %----------------------------------------------------------------------
    % Monomials
    N_2 = 2 * ceil(z_polesCorner_N/2);
    z_ast = 1e-12 + 0i*1;
    
    N_vec(j) = 2*N_1+2*N_2+1;
    
    %----------------------------------------------------------------------
    % System matrix
    A = GenSysMatrix_New(k,M,N_1,N_2,z_samp,z_samp_weights,z_poles,z_ast);
%     figure
%     imagesc(abs(A)); colorbar
    
%     A_ = A;
%     n_ = sqrt(sum(abs(A_).^2,1)); % Compute norms of columns
%     A_ = bsxfun(@rdivide,A_,n_); % Normalize M
%     n_ = reshape(n_,[],1); % Store column vector of norms
    
    %----------------------------------------------------------------------
    % Right hand side
    opt = {'ps',1};
    u_in_ = SelectIncidentField(k,opt);
    
    b = u_in_(z_samp);
    
    %----------------------------------------------------------------------
    % Solve system
    x = A\b;
%     x_reg = A_\b;
    
    %----------------------------------------------------------------------
    % Parse coefficients
    a_re = x(1:N_1);
    a_im = x(N_1+1:2*N_1);
    a = [a_re;a_im];
    b_re = x(2*N_1+1:2*N_1+1+N_2);
    b_im = x(2*N_1+1+N_2+1:end);
    
    %----------------------------------------------------------------------
    % Store L_inf boundary error, system matrix condition numbers, and
    % coefficient norms
    A_cond(j) = cond(A);
%     A__cond(j) = cond(A_);
    x_norm(j) = norm(x);
%     x_reg_norm(j) = norm(x_reg);
    err(j) = norm(A*x-b,2);
    fprintf(...
        'z_N: %d   ||Ax - b||_inf: %.2e   (cond(A): %.2e   sqrt(N): %.2f   MxN: %dx%d)\n',...
        z_polesCorner_N, err(j),cond(A),sqrt(N_vec(j)),M,N_vec(j));
end

%----------------------------------------------------------------------
% Plot convergence
%----------------------------------------------------------------------
sqrt_N = sqrt(N_vec);

hFig = figure;
set(hFig, 'Position', [100 100 1520 900])
subplot(3,2,1)
hold on, grid on
plot(sqrt_N,err,'b');
plot(sqrt_N,err,'b*');
set(gca,'yscale','log')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('$||Ax-b||_\infty$','Interpreter','latex')

subplot(3,2,3)
hold on, grid on
plot(sqrt_N,A_cond,'b');
plot(sqrt_N,A_cond,'b*');
% plot(sqrt_N,A__cond,'r');
% plot(sqrt_N,A__cond,'r*');
set(gca,'yscale','log')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('$cond(A)$','Interpreter','latex')

subplot(3,2,5)
hold on, grid on
plot(sqrt_N,x_norm,'b');
plot(sqrt_N,x_norm,'b*');
% plot(sqrt_N,x_reg_norm,'r');
% plot(sqrt_N,x_reg_norm,'r*');
set(gca,'yscale','log')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('$||x||_2$','Interpreter','latex')

%----------------------------------------------------------------------
% Evaluate on grid
%----------------------------------------------------------------------
x_Max = 1.0;
x_Min = -x_Max;

x_step = 0.020;
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
% hFig = figure;
% set(hFig, 'Position', [980 100 760 900])
subplot(3,2,2)
hold on, grid on
scatter(real(w_vec),imag(w_vec),'k*');
scatter(real(z_poles(:)),imag(z_poles(:)),'r.');
scatter(real(z_ast),imag(z_ast),'b*');
scatter(real(z_samp),imag(z_samp),'g+');
axis equal

subplot(3,2,4)
caxisVals = [];
GoodCAxis(real(u_tot));
PlotField(X,Y,real(u_tot),x_Max,caxisVals);

subplot(3,2,6)
caxisVals = [];
GoodCAxis(imag(u_tot));
PlotField(X,Y,imag(u_tot),x_Max,caxisVals);
