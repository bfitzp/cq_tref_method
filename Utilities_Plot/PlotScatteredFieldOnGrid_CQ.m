function PlotScatteredFieldOnGrid_CQ( ...
    ww,wr,wi,wc,scl, ...
    sig,z_0,c, ...
    p,time_N,T,dt,testFreqResults,...
    grid_N)

%----------------------------------------------------------------------
% Sample signal
%----------------------------------------------------------------------
t = linspace(0,T,time_N+1);
g = sig(t);

%----------------------------------------------------------------------
% Perform scaled FFT
%----------------------------------------------------------------------
R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));
h = bsxfun(@times,g,R.^(0:time_N));
h = fft(h,[],2);

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

%----------------------------------------------------------------------
% Generated observation grid
%----------------------------------------------------------------------
[x,y,z] = GenGrid(grid_N,wr,wi);

%----------------------------------------------------------------------
% Evaluate scattered field on observation grid
%----------------------------------------------------------------------
u_s = zeros(size(z(:),1),time_N+1);

for l=0:length(testFreqResults)-1
    k_l = 1i*p(R*omega^(-l))/(c*dt);
    
    c_adapt = testFreqResults{l+1}{1};
    nrmlz_adapt = testFreqResults{l+1}{2}.';
    pol_adapt = testFreqResults{l+1}{3};
    n_adapt = testFreqResults{l+1}{4};

    u_s(:,l+1) = fzeval_Helm_testFreq(nrmlz_adapt,k_l,z(:),wc,...
        c_adapt,pol_adapt,n_adapt);
    
    fprintf('l: %d   ---   k_l:   %.2e i%.2e   ---   norm(u_s): %.4e\n', ...
        l, real(k_l), imag(k_l), norm(u_s(:,l+1)));
end

u_s(:,time_N+2-(1:floor(time_N/2))) = conj(u_s(:,2:floor(time_N/2)+1));
u_s = real(ifft(u_s,[],2));
u_s = bsxfun(@times,u_s,R.^(-(0:time_N)));

%----------------------------------------------------------------------
% Evaluate incident field on observation grid
%----------------------------------------------------------------------
diffs = bsxfun(@minus,z(:),z_0);
dist = abs(diffs);

u_in = zeros(size(z(:),1),time_N+1);
for l=freqs_Range
    k_l = 1i*p(R*omega^(-l))/(c*dt);
    A = 1i/4*besselh(0,1,k_l*dist);
    u_in(:,l+1) = A*h(:,l+1);
    fprintf('l: %d   ---   k_l:   %.2e i%.2e   ---   norm(h): %.4e   ---   norm(A): %.4e   ---   norm(u_Bdy): %.4e\n', ...
        l, real(k_l), imag(k_l), norm(h(:,l+1)),norm(A),norm(u_in(:,l+1)));
end

u_in(:,time_N+2-(1:floor(time_N/2))) = conj(u_in(:,2:floor(time_N/2)+1));
u_in = ifft(u_in,[],2);
u_in = bsxfun(@times,u_in,R.^(-(0:time_N)));

% u_tot = u_in;
% u_tot = u_s___;
u_tot = u_s + u_in;

%----------------------------------------------------------------------
% Plot setup
%----------------------------------------------------------------------
LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
%     axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
% axes(PO,[.52 .34 .47 .56])

inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));

hFig = figure;
colorbar
% maxCLim = 1.1*max(real(u_tot(:,time_N)));
maxCLim = 0.08;
for t=1:time_N
    ut = u_tot(:,t);
    ut = real(reshape(ut,size(z)));
    
    ut(~inpolygonc(z,ww)) = nan;
     
    figure(hFig)
    clf
    imagesc(x,y,ut), colorbar, axis equal, hold on, set(gca,'ydir','normal')
    colormap(jet), caxis([-maxCLim,maxCLim])
    plot(ww,'-k',LW,1)
    % plot(pol,'.r',MS,6)
    set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6), axis(ax)
    hold on, grid on
    
    %     ylim([-maxYLim,maxYLim])
    %     ylim([-0.0574,0.0574])
    
    timeStr = sprintf('Time: %.2f', t*dt);
    title(timeStr);
    %         title(sprintf(strcat(titleStr, timeStr())));
    pause(0.01)
end


end

