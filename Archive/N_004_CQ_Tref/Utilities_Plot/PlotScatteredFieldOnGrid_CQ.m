function PlotScatteredFieldOnGrid_CQ( ...
        P,ww,wr,wi,wc,scl,n, ...
        sig,nrmlz_vec,z_0,c, ...
        coeffs,coeffs_fft,u_Bdy_coeffs,Z,pol, ...
        p,time_N,T,dt,k_threshold_Ind)

f = @(z,c,nrmlz,k) fzeval_Helm(nrmlz,k,z(:),wc,...               % vector and matrix inputs
    c,pol,n);
u_s_ = @(z,c,nrmlz,k) f(z,c,nrmlz,k);

inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

% Sample signal
t = linspace(0,T,time_N+1);
g = sig(t);
h = bsxfun(@times,g,R.^(0:time_N));    
h = fft(h,[],2);

if 1
    LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
    fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
    s_N = 100;
    sx = linspace(wr(1),wr(2),s_N); sy = linspace(wi(1),wi(2),s_N);
    [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
    ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
    axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
end

f = @(z,c,nrmlz,k) fzeval_Helm(nrmlz,k,z(:),wc,...
    c,pol,n);
u_s_ = @(z,c,nrmlz,k) f(z,c,nrmlz,k);

inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));

% Spatial information
diffs = bsxfun(@minus,zz(:),z_0);
dist = abs(diffs);


%----------------------------------------------------------------------
% Evaluate
%----------------------------------------------------------------------
u_in = zeros(size(zz(:),1),time_N+1);
u_s = zeros(size(zz(:),1),time_N+1);
% parfor l=freqs_Range
% parfor l=freqs_Range
fprintf('\n\n')
for l=freqs_Range
    if l > k_threshold_Ind
        continue;
    else
        k_l = 1i*p(R*omega^(-l))/(c*dt);        
        u_s(:,l+1) = u_s_(zz,coeffs_fft(:,l+1),nrmlz_vec(:,l+1),k_l);
        
        fprintf('l: %d   ---   k_l:   %.2e i%.2e   ---   norm(u_s): %.4e\n', ...
            l, real(k_l), imag(k_l), norm(u_s(:,l+1)));
    end
end

u_s(:,time_N+2-(1:floor(time_N/2))) = conj(u_s(:,2:floor(time_N/2)+1));
u_s = real(ifft(u_s,[],2));
u_s = bsxfun(@times,u_s,R.^(-(0:time_N)));


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
% 
% u_tot = u_in;
% u_tot = u_s;
u_tot = -u_s + u_in;


axes(PO,[.52 .34 .47 .56])
% levels = linspace(min(G),max(G),20);

figure
colorbar
% maxCLim = 1.1*max(u(:,time_N));
maxCLim = 0.1;
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


end

