%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Goals:
%   1. Fix issue with incident field generation - recover deltaBEM result
%   2. Tidy
%       Deal with fft signals throughout
%           Convert to time domain when needed
%   3. Robust adaptive error (Next step)
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
% time_N = 150;
% T = 1.25;
time_N = 75;                    % Number of time steps
T = 0.625;                      % End time

dt = T/time_N;

p = @(z) 1.5-2*z+0.5*z.^2;      % ODE solver/tranfer function argument

sigOpt = 3;                     % Select signal to generate
geomOpt = 'square' ;            % Select geometry to generate
% geomOpt = 'l_shape';            % Select geometry to generate
% geomOpt = 'l_shape_nc';         % Select geometry to generate

grid_N = 100;                   % DoFs for observation grid

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
z_0 = 0.2 + 0.0i;

Z_hr_N = 1000;

poles_N = 50;
[Z,Z_hr,pol,d,tt,pt,g, ...
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
[coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs, testFreqResults] = SolveCQ( ...
    P,u_in_bdy_fft,p,c, ...
    Z,pol,d, ...
    T,time_N,dt,k_threshold_Ind, ...
    tt,pt,g, ...
    ww,wc,nw, ...
    z_0, ...
    Z_hr, ...
    u_in_bdy_fft_hr);


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
        P,ww,wr,wi,wc,scl,n, ...
        sig,nrmlz_vec,z_0,c, ...
        coeffs_fft,u_Bdy_coeffs ,Z,pol, ...
        p,time_N,T,dt,k_threshold_Ind,testFreqResults, ...
        grid_N)
end





% %----------------------------------------------------------------------
% % Generalizing Trefethen's Method to the time domain
% %----------------------------------------------------------------------
% % Goals:
% %   1. Tidy
% %   2. Robust adaptive error (Next step)
% %----------------------------------------------------------------------
% clc
% clear all
% close all
% 
% %----------------------------------------------------------------------
% % Includes
% %----------------------------------------------------------------------
% addpath('Utilities_CQ')
% addpath('Utilities_General')
% addpath('Utilities_Geometry')
% addpath('Utilities_Plot')
% 
% %----------------------------------------------------------------------
% % Plotting parameters
% %----------------------------------------------------------------------
% 
% %----------------------------------------------------------------------
% % Signal
% %----------------------------------------------------------------------
% sigOpt = 2;
% sig = GenSignal(sigOpt);
% 
% %----------------------------------------------------------------------
% % Time discretization parameters
% %----------------------------------------------------------------------
% % time_N = 300;
% % time_N = 200;
% % time_N = 300;
% % T = 2.5;
% % time_N = 150;
% % T = 1.25;
% time_N = 75;
% T = 0.625;
% dt = T/time_N;
% 
% %----------------------------------------------------------------------
% % ODE solver/Tranfer function argument
% %----------------------------------------------------------------------
% p = @(z) 1.5-2*z+0.5*z.^2;
% 
% %----------------------------------------------------------------------
% % Plot the frequencies associated with the ODE solver
% %----------------------------------------------------------------------
% doPlotFreqs_CQ = 0;
% if doPlotFreqs_CQ
%     PlotCQFreqs(sig,T,time_N,p)
% end
% 
% %----------------------------------------------------------------------
% % Geometry
% %----------------------------------------------------------------------
% geomOpt = 'square';
% % geomOpt = 'l_shape';
% % geomOpt = 'l_shape_nc';
% P = SelectGeometry(geomOpt);
% P_N = length(P);
% 
% %----------------------------------------------------------------------
% % Generate sample points and poles
% %----------------------------------------------------------------------
% z_0 = 0.2 + 0.0i;
% 
% Z_hr_N = 1000;
% 
% poles_N = 50;
% [Z,Z_hr,pol,d,tt,pt,g, ...
%     ww,wr,wi,wc,nw,scl] = GenSamplePointsAndBdyData(...
%     z_0, P, poles_N, Z_hr_N);
% 
% sx = linspace(wr(1),wr(2),200); sy = linspace(wi(1),wi(2),200);
%     [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
% src = [real(z_0) imag(z_0)];
% % u_in = cylindricalwave(sig,1,src,[real(zz(:)) imag(zz(:))],0*[real(zz(:)) imag(zz(:))],T,time_N,0);
% [u_in, u_in_fft] = cylindricalwave(sig,1,src,[real(Z) imag(Z)],0*[real(Z) imag(Z)],T,time_N,0);
% 
% 
% %----------------------------------------------------------------------
% % Propagate incident field to boundary
% %----------------------------------------------------------------------
% c = 1;
% [u_inc_Bdy, u_inc_Bdy_fft, k_threshold_Ind] = PropagteIncidentField(sig,p,c,z_0,Z,T,time_N);
% 
% doPlot_u_inc_Bdy_CQ = 0;
% if doPlot_u_inc_Bdy_CQ
%     Plot_BdyField_CQ(P,Z,[],real(u_inc_Bdy),1,time_N,dt);
% end
% 
% % High resolution incident field
% % b_ = cell(nw);
% % b_fft_ = cell(nw);
% % for k = 1:nw
% %     [b_{k},b_fft_{k},~] = PropagteIncidentField(sig,p,c,z_0,Z_{k},T,time_N);
% % end
% [u_in_hr,u_in_fft_hr,~] = PropagteIncidentField(sig,p,c,z_0,Z_hr,T,time_N);
% 
% 
% 
% % if doPlot_u_inc_Bdy_CQ
% %     Plot_BdyField_CQ(P,Z,[],real(u_in),1,time_N,dt);
% % end
% 
% %----------------------------------------------------------------------
% % Solve time-domain boundary integral equation
% %----------------------------------------------------------------------
% % [coeffs, coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs] = SolveCQ( ...
% %     u_inc_Bdy,u_inc_Bdy_fft,p,c, ...
% %     Z,pol,d, ...
% %     T,time_N,dt,k_threshold_Ind, ...
% %     tt,pt,g,...
% %     ww,wc,nw);
% 
% [coeffs, coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs, testFreqResults] = SolveCQ( ...
%     P,u_in,u_in_fft,p,c, ...
%     Z,pol,d, ...
%     T,time_N,dt,k_threshold_Ind, ...
%     tt,pt,g, ...
%     ww,wc,nw, ...
%     z_0, ...
%     Z_hr, ...
%     u_in_hr,u_in_fft_hr);
% 
% 
% % doPlotScatteredFieldOnBdy_CQ = 0;
% % if doPlotScatteredFieldOnBdy_CQ
% %     PlotScatteredFieldOnBdy_CQ(...
% %         P,ww,wr,wi,wc,scl,n, ...
% %         nrmlz_vec, ...
% %         coeffs,coeffs_fft,u_Bdy_coeffs,Z,pol, ...
% %         p,time_N,T,dt,k_threshold_Ind)
% % end
% 
% %----------------------------------------------------------------------
% % Plot field on grid
% %----------------------------------------------------------------------
% doPlotScatteredFieldOnGrid_CQ = 1;
% if doPlotScatteredFieldOnGrid_CQ
%     PlotScatteredFieldOnGrid_CQ( ...
%         P,ww,wr,wi,wc,scl,n, ...
%         sig,nrmlz_vec,z_0,c, ...
%         coeffs,coeffs_fft,u_Bdy_coeffs ,Z,pol, ...
%         p,time_N,T,dt,k_threshold_Ind,testFreqResults)
% end