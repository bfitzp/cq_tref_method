%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
clc
clear all
close all

%----------------------------------------------------------------------
% Includes
%----------------------------------------------------------------------
AddPaths();

%----------------------------------------------------------------------
% Simulation parameters
%----------------------------------------------------------------------
c = 1;                          % Speed of sound

% time_N = 300;
% T = 2.5;
% time_N = 240;
% T = 2.00;
time_N = 150;
T = 1.25;
% time_N = 75;                    % Number of time steps
% T = 0.625;                      % End time

dt = T/time_N;

p = @(z) 1.5-2*z+0.5*z.^2;      % ODE solver/tranfer function argument

sigOpt = 3;                     % Select signal to generate
% geomOpt = 'square' ;            % Select geometry to generate
geomOpt = 'l_shape';            % Select geometry to generate
% geomOpt = 'l_shape_nc';         % Select geometry to generate

grid_N = 150;                   % DoFs for observation grid

errTol = 1e-2;                 % Desired error tolerance

%----------------------------------------------------------------------
% Generate signal and geometry
%----------------------------------------------------------------------
sig = GenSignal(sigOpt);
P = GenGeometry(geomOpt);

%----------------------------------------------------------------------
% Plot the frequencies associated with the ODE solver
%----------------------------------------------------------------------
doPlotFreqs_CQ = 0;
if doPlotFreqs_CQ
    PlotCQFreqs(sig,T,time_N,p)
end

%----------------------------------------------------------------------
% Generate sample points and poles
%----------------------------------------------------------------------
z_0 = 0.2 + 0.0i; % Trouble
% z_0 = 0.2 + 0.3i;

Z_hr_N = 1000;

poles_N = 50;
[Z,Z_hr, ...
    ww,wr,wi,wc,nw,scl] = GenSamplePointsAndBdyData(...
    z_0, P, poles_N, Z_hr_N);

%----------------------------------------------------------------------
% Propagate incident field to boundary
%----------------------------------------------------------------------
[u_in_bdy_fft, k_threshold_Ind] = PropagteIncidentField( ...
    sig,p,c,z_0,Z,T,time_N,dt);

[u_in_bdy_fft_hr, ~] = PropagteIncidentField( ...
    sig,p,c,z_0,Z_hr,T,time_N,dt);

doPlot_u_inc_Bdy_CQ = 0;
if doPlot_u_inc_Bdy_CQ
    Plot_BdyField_CQ(P,Z,u_in_bdy_fft,time_N,dt);
end

%----------------------------------------------------------------------
% Solve time-domain boundary integral equation
%----------------------------------------------------------------------
testFreqResults = SolveCQ( ...
    P, ...
    u_in_bdy_fft_hr, ...
    p,c, ...
    time_N,dt, ...
    Z_hr, ...
    errTol);

% doPlotScatteredFieldOnBdy_CQ = 0;
% if doPlotScatteredFieldOnBdy_CQ
%     PlotScatteredFieldOnBdy_CQ(...
%         P,ww,wr,wi,wc,scl,n, ...
%         nrmlz_vec, ...
%         coeffs,coeffs_fft,u_Bdy_coeffs,Z,pol, ...
%         p,time_N,T,dt,k_threshold_Ind)
% end

%----------------------------------------------------------------------
% Plot field on grid
%----------------------------------------------------------------------
doPlotScatteredFieldOnGrid_CQ = 1;
if doPlotScatteredFieldOnGrid_CQ
    PlotScatteredFieldOnGrid_CQ( ...
        ww,wr,wi,wc,scl, ...
        sig,z_0,c, ...
        p,time_N,T,dt,testFreqResults, ...
        grid_N)
end