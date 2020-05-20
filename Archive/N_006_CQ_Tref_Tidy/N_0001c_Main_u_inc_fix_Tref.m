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
% sig = @(t) sin(2*t).^5.*(t>=0);
sig = @(t) sin(16*t).^5.*(t>=0);
% sig = @(t) sin(32*t).^5.*(t>=0);

x_0 = 0.5-0.5i;

% --------------------------------------------------
% Time discretization parameters
% --------------------------------------------------
% time_N = 300;
% T = 2.5;
% time_N = 150;
% T = 1.25;
% time_N = 200;
% T = 2.00;
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
% geomOpt = 'square';
geomOpt = 'l_shape';
P = SelectGeometry(geomOpt);
P_N = length(P);

%----------------------------------------------------------------------
% Generate sample points and poles
%----------------------------------------------------------------------
z_0 = 0.5 + 0.0i;
% z_0 = 0.8 + 0.3i;

poles_N = 60;
[Z,pol,d,tt,pt,g, ...
    ww,wr,wi,wc,nw,scl] = GenSamplePointsAndBdyData(...
    z_0, P, poles_N);

c = 1;





LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
sx = linspace(wr(1),wr(2),200); sy = linspace(wi(1),wi(2),200);
[xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';

inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));

% Spatial information
diffs = bsxfun(@minus,zz(:),z_0);
dist = abs(diffs);

% Sample signal
t = linspace(0,T,time_N+1);
g = sig(t);
dt = T/time_N;
R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

h = bsxfun(@times,g,R.^(0:time_N));    
h = fft(h,[],2);



u_in = zeros(size(zz(:),1),time_N+1);
for l=freqs_Range
    k_l = 1i*p(R*omega^(-l))/(c*dt);

    A = 1i/4*besselh(0,1,k_l*dist);
    u_in(:,l+1) = A*h(:,l+1); 
end


u_in(:,time_N+2-(1:floor(time_N/2))) = conj(u_in(:,2:floor(time_N/2)+1));

u_in = ifft(u_in,[],2);
u_in = bsxfun(@times,u_in,R.^(-(0:time_N)));

u_tot = u_in;



axes(PO,[.52 .34 .47 .56])
% levels = linspace(min(G),max(G),20);

figure
colorbar
% maxCLim = 1.1*max(u(:,time_N));
maxCLim = 0.2;
for t=1:time_N
    ut = u_tot(:,t);
    ut = real(reshape(ut,size(zz)));

    ut(~inpolygonc(zz,ww)) = nan;

    clf
    imagesc(sx,sy,ut), colorbar, axis equal, hold on, set(gca,'ydir','normal')
    colormap(jet), caxis([-maxCLim,maxCLim])
    plot(ww,'-k',LW,1), plot(pol,'.r',MS,6)
    set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6), axis(ax)
    hold on, grid on
    
%     ylim([-maxYLim,maxYLim])
%     ylim([-0.0574,0.0574])
    
    timeStr = sprintf('Time: %.2f', t*dt);
    title(timeStr);
    %         title(sprintf(strcat(titleStr, timeStr())));
    pause(0.01)
end
