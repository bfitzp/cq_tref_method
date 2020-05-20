function PlotScatteredFieldOnBdy_CQ(...
    P,ww,wr,wi,wc,scl,n, ...
    nrmlz_vec, ...
    coeffs,coeffs_fft,u_Bdy_coeffs,Z,pol, ...
    p,time_N,dt,k_threshold_Ind)

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

if 1
    LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
    fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
    sx = linspace(wr(1),wr(2),200); sy = linspace(wi(1),wi(2),200);
    [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
    ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
    axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
end

% f = @(z,c,nrmlz,k) reshape(fzeval_Helm(nrmlz,k,z(:),wc,...               % vector and matrix inputs
%     c,pol,n),size(z));
f = @(z,c,nrmlz,k) fzeval_Helm(nrmlz,k,z(:),wc,...               % vector and matrix inputs
    c,pol,n);
u_ = @(z,c,nrmlz,k) f(z,c,nrmlz,k);

inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));  

%----------------------------------------------------------------------
% Evaluate
%----------------------------------------------------------------------
u = zeros(size(Z,1),time_N+1);
% parfor l=freqs_Range
for l=freqs_Range
    if l > k_threshold_Ind/1
        continue;
    else
        fprintf('Evaluating solution: Step %d of %d\n',l,k_threshold_Ind);
        
        k_l = 1i*p(R*omega^(-l))/dt;
        
        u(:,l+1) = u_(Z,coeffs_fft(:,l+1),nrmlz_vec(:,l+1),k_l);

%         u(:,l+1) = uu;
    end
end

u(:,time_N+2-(1:floor(time_N/2))) = conj(u(:,2:floor(time_N/2)+1));

u = real(ifft(u,[],2));
u = bsxfun(@times,u,R.^(-(0:time_N)));

Plot_BdyField_CQ(P,Z,[],real(u),1,time_N,dt);

end