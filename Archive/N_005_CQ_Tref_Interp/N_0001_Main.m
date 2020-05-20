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
close all

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
sigOpt = 2;
switch(sigOpt)
    case 1
        sig = @(t) sin(2*t).^5.*(0<=t);
    case 2
        sig = @(t) sin(16*t).^5.*(0<=t).*(t<=0.40);
    case 3
        sig = @(t) sin(32*t).^5.*(0<=t);
    case 4
        sig = @(t) sin(32*t).^5.*(0<=t).*(t<=0.49);
end

% figure
% plot(0:0.001:1.25,sig(0:0.001:1.25));

% x_0 = 0.5-0.5i;

% --------------------------------------------------
% Time discretization parameters
% --------------------------------------------------
% time_N = 300;
% time_N = 200;
% time_N = 300;
% T = 2.5;
% time_N = 150;
% T = 1.25;
time_N = 75;
T = 0.625;
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
% geomOpt = 'l_shape';
% geomOpt = 'l_shape_nc';
P = SelectGeometry(geomOpt);
P_N = length(P);

%----------------------------------------------------------------------
% Generate sample points and poles
%----------------------------------------------------------------------
z_0 = 0.2 + 0.0i;
% z_0 = 0.8 + 0.3i;

Z_hr_N = 1000;

poles_N = 50;
[Z,Z_hr,pol,d,tt,pt,g, ...
    ww,wr,wi,wc,nw,scl] = GenSamplePointsAndBdyData(...
    z_0, P, poles_N, Z_hr_N);

sx = linspace(wr(1),wr(2),200); sy = linspace(wi(1),wi(2),200);
    [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
src = [real(z_0) imag(z_0)];
% u_in = cylindricalwave(sig,1,src,[real(zz(:)) imag(zz(:))],0*[real(zz(:)) imag(zz(:))],T,time_N,0);
[u_in, u_in_fft] = cylindricalwave(sig,1,src,[real(Z) imag(Z)],0*[real(Z) imag(Z)],T,time_N,0);


%----------------------------------------------------------------------
% Propagate incident field to boundary
%----------------------------------------------------------------------
c = 1;
[u_inc_Bdy, u_inc_Bdy_fft, k_threshold_Ind] = PropagteIncidentField(sig,p,c,z_0,Z,T,time_N);

doPlot_u_inc_Bdy_CQ = 0;
if doPlot_u_inc_Bdy_CQ
    Plot_BdyField_CQ(P,Z,[],real(u_inc_Bdy),1,time_N,dt);
end

% High resolution incident field
% b_ = cell(nw);
% b_fft_ = cell(nw);
% for k = 1:nw
%     [b_{k},b_fft_{k},~] = PropagteIncidentField(sig,p,c,z_0,Z_{k},T,time_N);
% end
[u_in_hr,u_in_fft_hr,~] = PropagteIncidentField(sig,p,c,z_0,Z_hr,T,time_N);



% if doPlot_u_inc_Bdy_CQ
%     Plot_BdyField_CQ(P,Z,[],real(u_in),1,time_N,dt);
% end

%----------------------------------------------------------------------
% Solve time-domain boundary integral equation
%----------------------------------------------------------------------
% [coeffs, coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs] = SolveCQ( ...
%     u_inc_Bdy,u_inc_Bdy_fft,p,c, ...
%     Z,pol,d, ...
%     T,time_N,dt,k_threshold_Ind, ...
%     tt,pt,g,...
%     ww,wc,nw);

[coeffs, coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs, testFreqResults] = SolveCQ( ...
    P,u_in,u_in_fft,p,c, ...
    Z,pol,d, ...
    T,time_N,dt,k_threshold_Ind, ...
    tt,pt,g, ...
    ww,wc,nw, ...
    z_0, ...
    Z_hr, ...
    u_in_hr,u_in_fft_hr);


% doPlotScatteredFieldOnBdy_CQ = 0;
% if doPlotScatteredFieldOnBdy_CQ
%     PlotScatteredFieldOnBdy_CQ(...
%         P,ww,wr,wi,wc,scl,n, ...
%         nrmlz_vec, ...
%         coeffs,coeffs_fft,u_Bdy_coeffs,Z,pol, ...
%         p,time_N,T,dt,k_threshold_Ind)
% end


%----------------------------------------------------------------------
% Plot geometry
%----------------------------------------------------------------------

%----------------------------------------------------------------------
% Plot field on grid
%----------------------------------------------------------------------
% n = 120;
% nrmlz_vec = [];
doPlotScatteredFieldOnGrid_CQ = 1;
if doPlotScatteredFieldOnGrid_CQ
    PlotScatteredFieldOnGrid_CQ( ...
        P,ww,wr,wi,wc,scl,n, ...
        sig,nrmlz_vec,z_0,c, ...
        coeffs,coeffs_fft,u_Bdy_coeffs ,Z,pol, ...
        p,time_N,T,dt,k_threshold_Ind,testFreqResults)
%     PlotScatteredFieldOnGrid_CQ( ...
%         P,ww,wr,wi,wc,scl,n, ...
%         sig,nrmlz_vec,z_0,c, ...
%         [],[],[],Z,pol, ...
%         p,time_N,T,dt,k_threshold_Ind)
end