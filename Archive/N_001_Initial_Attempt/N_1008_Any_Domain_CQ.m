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
addpath('Utilities_CQ')
addpath('Utilities_Plot')

%----------------------------------------------------------------------
% Plotting parameters
%----------------------------------------------------------------------
dx = 0.02;
% cvg_FigPos = [1920 100 760 900];
geom_FigPos = [0 100 760 900];

% --------------------------------------------------
% Input signal
% --------------------------------------------------
signal = @(t) sin(16*t).^5.*(t>=0);

x_0 = 0.5-0.5i;

% --------------------------------------------------
% Time discretization parameters
% --------------------------------------------------
time_N = 200;
T = 2;
dt = T/time_N;

% --------------------------------------------------
% ODE solver/Tranfer function argument
% --------------------------------------------------
% p_trap = @(z) 2*(1-z)/(1+z); % Also try BDF2 - Does deltaBEM use both?
p_bdf2 = @(z) 1.5-2*z+0.5*z.^2;
p = p_bdf2;

% --------------------------------------------------
% Plot the frequencies associated with the ODE solver
% --------------------------------------------------
doPlotFreqs_CQ = 1;
if doPlotFreqs_CQ
%     PlotCQFreqs(signal,T,M,{p_bdf2,p_trap})
    PlotCQFreqs(signal,T,time_N,{p_bdf2})
end

%----------------------------------------------------------------------
% Corners
%----------------------------------------------------------------------
geomOpt = 4;
[w, w_bndBox, cornerInds] = SelectGeometry(geomOpt);

%----------------------------------------------------------------------
% Poles
%----------------------------------------------------------------------
sigma = 0.7;
z_polesCorner_N = 120;
z_polesOffset = 0.00;

[z_poles,N_1]...
    = GenPoles_New(...
    z_polesCorner_N,z_polesOffset,sigma,...
    w, cornerInds);

%----------------------------------------------------------------------
% Monomials
N_2 = 1 * ceil(z_polesCorner_N/2 * 2);
z_ast = 0.0 + 0.0i;

%----------------------------------------------------------------------
% Sample points
%----------------------------------------------------------------------
edgePanels_N = 480;
refine_N = 1;
[z_samp,z_samp_weights,cornerInds,M] = GenSamplePoints_New(w,edgePanels_N,refine_N);

%----------------------------------------------------------------------
% Plot geometry
%----------------------------------------------------------------------
PlotGeometry(geom_FigPos,w,[],z_samp,0);

% --------------------------------------------------
% Compute the frequency range we need to work with
%   after discarding frequencies with negligible boundary norms
% Actually we should be doing this for Fourier modes?
% --------------------------------------------------
incFieldOpt = {'ps',1};
[~, z_0] = SelectIncidentField(1,incFieldOpt);
threshold = 1e-15;
[k_threshold,k_threshold_Ind,k_re_max] = ...
    ComputeFreqThreshold(threshold,time_N,dt,p,z_samp,incFieldOpt,0);

fprintf('k_threshold: %.4f %.4fi   ...   k_threshold_Ind: %d\n',...
    abs(real(k_threshold)), imag(k_threshold), k_threshold_Ind);

% --------------------------------------------------
% Propagate incident field to boundary
% --------------------------------------------------
c = 1;
[u_inc_Bdy, u_inc_Bdy_fft]  = PropagteCylindricalWave(signal,p,c,z_0,z_samp,T,time_N,k_threshold_Ind);
% u_inc_bdy_fast = PropagteCylindricalWave_fast(signal,p,c,z_0,z_samp,T,time_N,k_threshold_Ind);

doPlot_u_inc_Bdy_CQ = 1;
if doPlot_u_inc_Bdy_CQ
    Plot_BdyField_CQ(z_samp,cornerInds,real(u_inc_Bdy),1,time_N,dt);
%     Plot_BdyField_CQ(z_samp,cornerInds,real(u_inc_bdy_fast),1,time_N,dt);
end

% --------------------------------------------------
% Solve time-domain boundary integral equation
% --------------------------------------------------
% coeffs_CQ = SolveCQ(M_modes_Max,bdyData_r_theta,a,O,kappa,p);
[coeffs, coeffs_fft, AScale] = SolveCQ(u_inc_Bdy,p,c,...
    z_samp,z_samp_weights,z_poles,z_ast,N_2,...
    T,time_N,dt,k_threshold_Ind);

doPlotScatteredFieldOnBdy_CQ = 1;
if doPlotScatteredFieldOnBdy_CQ
    PlotScatteredFieldOnBdy_CQ(...
    coeffs,coeffs_fft,w,z_samp,z_poles,z_ast,...
    cornerInds,N_2,...
    p,time_N,dt,k_threshold_Ind,AScale)
end

% PlotFieldOnGrid_CQ(coeffs,p,c,...
%     w,w_bndBox,z_poles,z_ast,N_2,dx,...
%     T,time_N,dt,k_threshold_Ind)

% --------------------------------------------------
% % Compare negative incident field and scattered field on boundary
% % --------------------------------------------------
% doPlot_uBdy_CQ = 1;
% if doPlot_uBdy_CQ
%     u_s_bdy = Eval_Field_CQ(coeffs_CQ);
%     Plot_BdyField_CQ(z_samp,real(u_s_bdy),time_N,dt);
%     Plot_BdyField_CQ(z_samp,real(-u_inc_bdy),time_N,dt);
% end

%----------------------------------------------------------------------
% Plot geometry
%----------------------------------------------------------------------
PlotGeometry(geom_FigPos,w,z_poles,z_samp,z_ast)

%----------------------------------------------------------------------
% Plot field on grid
%----------------------------------------------------------------------
PlotFieldOnGrid_CQ(coeffs,p,c,...
    w,w_bndBox,z_poles,z_ast,N_2,dx,...
    T,time_N,dt,k_threshold_Ind)
