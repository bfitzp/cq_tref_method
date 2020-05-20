%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Goals:
%   1. Sampling that is independent of edge length but dependent on
%       wavenumber
%   2. Propagate incident wave to the boundary
%----------------------------------------------------------------------
clc
clear all
% close all

%----------------------------------------------------------------------
% Includes
%----------------------------------------------------------------------
addpath('Utilities_CQ')
addpath('Utilities_General')
addpath('Utilities_Geometry')
addpath('Utilities_Plot')

%----------------------------------------------------------------------
% Plotting parameters
%----------------------------------------------------------------------

% --------------------------------------------------
% Input signal
% --------------------------------------------------
sig = @(t) sin(16*t).^5.*(t>=0);

x_0 = 0.5-0.5i;

% --------------------------------------------------
% Time discretization parameters
% --------------------------------------------------
% time_N = 200;
time_N = 100;
T = 2;
dt = T/time_N;

% --------------------------------------------------
% ODE solver/Tranfer function argument
% --------------------------------------------------
p = @(z) 1.5-2*z+0.5*z.^2;

% --------------------------------------------------
% Plot the frequencies associated with the ODE solver
% --------------------------------------------------
doPlotFreqs_CQ = 0;
if doPlotFreqs_CQ
    PlotCQFreqs(sig,T,time_N,p)
end

% --------------------------------------------------
% Get CQ frequency range and Nyquist frequency
% --------------------------------------------------
[k_Range, k_Nyq] = GetCQFreqRange(time_N,dt,p);

%----------------------------------------------------------------------
% Geometry
%----------------------------------------------------------------------
geomOpt = 'square';
P = SelectGeometry(geomOpt);

vararg = {};
[g, w, ww, pt, dw, tol, steps, plots, ...
    slow, rel, arnoldi, aaaflag] = ...
    parseinputs_(P,k_Nyq,vararg{:});
Zplot = ww;
nw = length(w);
wr = sort(real(ww)); wr = wr([1 end]);
wi = sort(imag(ww)); wi = wi([1 end]);
wc = mean(wr+1i*wi);
scl = max([diff(wr),diff(wi)]);


for k = 1:nw
    forward = pt{k}(.01*dw(k)) - w(k);            % small step toward next corner
    j = mod(k-2,nw)+1;
    backward = pt{j}(.99*dw(j)) - w(k);           % small step toward last corner
    tmp = 1i*backward*sqrt(-forward/backward);
    outward(k) = tmp/abs(tmp);                    % outward direction from corner
end

%----------------------------------------------------------------------
% Generate sample points and boundary data
%----------------------------------------------------------------------
polmax = 100;
[Z,G,tt] = GenSamplePointsAndBdyData(...
    w,wc,wr,wi,ww,dw,nw,pt,outward,g,polmax,scl,k_Nyq);

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

% --------------------------------------------------
% Propagate incident field to boundary
% --------------------------------------------------
c = 1;
[u_inc_Bdy, u_inc_Bdy_fft]  = PropagteCylindricalWave(sig,p,c,z_0,Z,T,time_N,k_threshold_Ind);

doPlot_u_inc_Bdy_CQ = 1;
if doPlot_u_inc_Bdy_CQ
    Plot_BdyField_CQ(z_samp,cornerInds,real(u_inc_Bdy),1,time_N,dt);
end

return









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