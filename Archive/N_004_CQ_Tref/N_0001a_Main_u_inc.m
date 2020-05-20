%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Goals:
%   1. Remove convergence process from Trefethen's method and
%       and use fixed DoF
%   2. CQ enable the result
%   3. Plot u_inc on boundary
%   3. Plot u_s on boundary
%   4. Plot u_s on grid
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
time_N = 150;
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

%----------------------------------------------------------------------
% Geometry
%----------------------------------------------------------------------
geomOpt = 'square';
P = SelectGeometry(geomOpt);
P_N = length(P);

%----------------------------------------------------------------------
% Generate sample points and poles
%----------------------------------------------------------------------
z_0 = 0.5 - 0.0i;
poles_N = 40;
[Z,pol,d,tt,pt,g, ...
    ww,wr,wi,wc,nw,scl] = GenSamplePointsAndBdyData(...
    z_0, P, poles_N);
%----------------------------------------------------------------------
% Propagate incident field to boundary
%----------------------------------------------------------------------
c = 1;
[u_inc_Bdy, u_inc_Bdy_fft, k_threshold_Ind] = PropagteIncidentField_u_inc(sig,p,c,z_0,Z,T,time_N);

doPlot_u_inc_Bdy_CQ = 1;
if doPlot_u_inc_Bdy_CQ
    Plot_BdyField_CQ(P,Z,[],real(u_inc_Bdy),1,time_N,dt);
end

