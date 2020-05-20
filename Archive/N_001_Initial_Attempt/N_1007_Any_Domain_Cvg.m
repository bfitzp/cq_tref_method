%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Convergence of error for frequency domain with square geometry
% as number of degrees of freedom increase (poles and monomials)
%----------------------------------------------------------------------
clc
clear all
close all

%----------------------------------------------------------------------
% Includes
%----------------------------------------------------------------------
addpath('Utilities')
addpath('Utilities_Plot')

%----------------------------------------------------------------------
% Plotting parameters
%----------------------------------------------------------------------
dx = 0.01;
cvg_FigPos = [100 100 760 900];
geomFields_FigPos = [900 100 760 900];

%----------------------------------------------------------------------
% Corners
%----------------------------------------------------------------------
opts = 1;
[w, w_bndBox, cornerInds] = SelectGeometry(opts);

%----------------------------------------------------------------------
% Wavenumber
%----------------------------------------------------------------------
k = 2 + 0i;

%----------------------------------------------------------------------
% Poles
%----------------------------------------------------------------------
z_polesCorner_Min = 15;
z_polesCorner_Max = 15;
z_polesCorner_Step = 4;
z_polesCorner_Range = ...
    z_polesCorner_Min:z_polesCorner_Step:z_polesCorner_Max;
z_polesCorner_Range_N = length(z_polesCorner_Range);

%----------------------------------------------------------------------
% Iterate with an increasing number of DoFs
%----------------------------------------------------------------------
A_cond = zeros(z_polesCorner_Range_N,1);
x_norm = zeros(z_polesCorner_Range_N,1);
bdyErr = zeros(z_polesCorner_Range_N,1);
N_vec = zeros(z_polesCorner_Range_N,1);

%----------------------------------------------------------------------
% Sample points
edgePanels_N = 80;
refine_N = 2;
[z_samp,z_samp_weights,cornerInds,M] = GenSamplePoints_New(w,edgePanels_N,refine_N);  
    
for j=1:z_polesCorner_Range_N
    %----------------------------------------------------------------------
    % Poles
    sigma = 0.7;
    z_polesCorner_N = z_polesCorner_Range(j);
    z_polesOffset = 0.00;
    
    [z_poles,N_1]...
        = GenPoles_New(...
        z_polesCorner_N,z_polesOffset,sigma,...
        w, cornerInds);
    
    %----------------------------------------------------------------------
    % Monomials
    N_2 = 1 * ceil(z_polesCorner_N/2 * 2);
    z_ast = 0.25 - 0.25i;
    
    N_vec(j) = 2*N_1+2*N_2+1;
    
    %----------------------------------------------------------------------
    % System matrix
    A = GenSysMatrix_MPS(k,N_2,z_samp,z_samp_weights,z_poles,z_ast);
    
    %----------------------------------------------------------------------
    % Right hand side
%     opt = {'ps',1};
    opt = {'fn',1};
    u_in_ = SelectIncidentField(k,opt);
    
    b = u_in_(z_samp);
%     disp(norm(b))
    
    %----------------------------------------------------------------------
    % Solve system
    coeffs = A\b;
%     coeffs = lsqminnorm(A,b);
%     coeffs = pinv(A)*b;
    
%     A_t_A = A.'*A;
%     A_t_b = A.'*b;
%     A_t_A = A'*A;
%     A_t_b = A'*b;
%     I = eye(size(A_t_A));
%     lam = 0.1; % 10e-3*trace(A_t_A)/size(A_t_A,2);
%     coeffs = (A_t_A + lam*I)\A_t_b;
    
    %---------------------------------------------------------------------
    % Store L_inf boundary error, system matrix condition numbers, and
    % coefficient norms
    [row, col] = find(isnan(A));
    A_cond(j) = cond(A);
    x_norm(j) = norm(coeffs);
    bdyErr(j) = norm(A*coeffs-b,2);
    fprintf(...
        'z_N: %d   ||Ax - b||_inf: %.2e   cond(A): %.2e   sqrt(N): %.2f   MxN: %dx%d\n',...
        z_polesCorner_N, bdyErr(j), A_cond(j), sqrt(N_vec(j)), M, N_vec(j));
    
%     %----------------------------------------------------------------------
%     % Plot geometry and fields
%     Plot_Geom_and_Fields(...
%         geomFields_FigPos,dx,coeffs,w,w_bndBox,...
%         k,z_poles,z_samp,z_samp_weights.^0,z_ast,N_2)
end

%----------------------------------------------------------------------
% Plot convergence
%----------------------------------------------------------------------
PlotConvergence(cvg_FigPos,N_vec,bdyErr,A_cond,x_norm)

%----------------------------------------------------------------------
% Plot geometry and fields
%----------------------------------------------------------------------
Plot_Geom_and_Fields(...
    geomFields_FigPos,dx,coeffs,w,w_bndBox,...
    k,z_poles,z_samp,z_samp_weights.^0,z_ast,N_2)
