% function [u, maxerr, f, Z, Zplot, A] = N_0003_Helmholtz_Eq_CQ_Export_(k_wave, P, z_0, varargin)
function [c, nrmlz, Z, pol, n, d] = N_0003_Helmholtz_Eq_CQ_Export_(...
    k_wave, P, z_0, Z_hr, b_hr, varargin)

doPlot = 0;
%----------------------------------------------------------------------
% Corners
%----------------------------------------------------------------------
% opts = 6;
% [w, w_bndBox, cornerInds] = SelectGeometry(opts);


%% Set up the problem
[g, w, ww, pt, dw, tol, steps, plots, ...        % parse inputs
    slow, rel, arnoldi, aaaflag] = ...
    parseinputs_(P,varargin{:});
Zplot = ww;
nw = length(w);                                  % number of corners
wr = sort(real(ww)); wr = wr([1 end]);
wi = sort(imag(ww)); wi = wi([1 end]);
wc = mean(wr+1i*wi);                             
scl = max([diff(wr),diff(wi)]);                 % for scale- and transl-invariance

%----------------------------------------------------------------------
% Modifications
%----------------------------------------------------------------------

tol = 1e-10;
tt_N = 40;

q = .5;                 % sets which corners get more poles
if slow == 1
    q = 0;
end                 % sets which corners get more poles
inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));            % complex variant of "inpolygon"

for k = 1:nw
    forward = pt{k}(.01*dw(k)) - w(k);            % small step toward next corner
    j = mod(k-2,nw)+1;
    backward = pt{j}(.99*dw(j)) - w(k);           % small step toward last corner
    tmp = 1i*backward*sqrt(-forward/backward);
    outward(k) = tmp/abs(tmp);                    % outward direction from corner
end
warn = warning('off','MATLAB:rankDeficientMatrix');  % matrices are ill-conditioned

%% Set up for plots
if plots
    LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
    fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
    sx = linspace(wr(1),wr(2),200); sy = linspace(wi(1),wi(2),200);
    [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
    ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
    axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
end

%% Main loop: increase number of poles until convergence ==============================
Nvec = []; errvec = []; condAvec = []; coeffvec = [];
tic
errk = ones(nw,1);                               % max error near each corner
nkv = zeros(nw,1);                               % no. of poles at each corner
minstepno = 1;
maxstepno = 35; err0 = Inf;

hFig_G = figure(1);
hFig_geom = figure(2);

firstFlag = 1;
for stepno = minstepno:maxstepno
    % Fix poles and sample pts on bndry.  Side k means side from corner k to k+1.
    Z = [];           % col vector of sample points on boundary
    G = [];           % col vector of boundary values at these points
    T = [];           % col vector of unit tangent vectors at these points
    pol = [];         % row vector of poles of the rational approximation
    J = [];           % row vector of indices of which corner each pole belongs to
    d = [];           % row vector of distances from poles to their corners
    tt = cell(nw,1);  % cell array of distances of sample points along each side
    
    %----------------------------------------------------------------------
    % Generate poles
    %----------------------------------------------------------------------
    for k = 1:nw
        nk = nkv(k);                                  % no. of poles at this corner
        sk = sqrt(1:nk) - sqrt(nk);
        dk = exp(2*sk); dk = scl*dk;                  % stronger clustering near corner
        dk = dk(dk>1e-15*scl);                        % remove poles too close to corner
        polk = w(k) + outward(k)*dk;                  % poles near this corner
        ii = find(inpolygonc(polk(dk>1e-12*scl),ww),1); % work around inaccuracy
        if length(ii)>0                               % don't allow poles in Omega
            dk = dk(1:ii-2); polk = polk(1:ii-2);
        end
        pol = [pol polk]; d = [d dk];
        
        J = [J k*ones(1,length(dk))];
        dvec = [(1/3)*dk (2/3)*dk dk];                % finer pts for bndry sampling
        tt{k} = [tt{k} dvec(dvec<dw(k)) ...           % add clustered pts near corner
            linspace(0,dw(k),max(tt_N,nk))];            % additional pts along side
        j = mod(k-2,nw)+1;                            % index of last corner
        tt{j} = [tt{j} dw(j)-dvec(dvec<dw(j))];       % likewise in other direction
    end
    
    figure(hFig_geom);
    clf
    hold on, grid on
    plot(ww,'-k',LW,1), plot(pol,'*r',MS,6), plot(Z,'+m',MS,6),
    set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6)
    axis(ax), axis equal
%     xlim([-0.1, 0.1]), ylim([-0.1, 0.1])
    xlim([-3, 3]), ylim([-3, 3])

    %----------------------------------------------------------------------
    % Generate sample points and boundary data
    %----------------------------------------------------------------------
    for k = 1:nw
        tt{k} = sort(tt{k}(:));
        tk = tt{k}; pk = pt{k};                       % abbrevations
        Z = [Z; pk(tk)];                              % sample pts on side k
        h = 1e-4;                                     % 4-pt trapezoidal rule
    end
    
    Z_hr_N = 1000;
    for k = 1:nw
        tt{k} = sort(tt{k}(:));
        tk = tt{k}; pk = pt{k};  
        pts = pk(tk);
        ind_0 = (k-1)*Z_hr_N+1;
        ind_1 = k*Z_hr_N;
        G = [G; InterpolateIncField(pts,Z_hr(ind_0:ind_1),b_hr(ind_0:ind_1))];
    end
                                % normalize tangent vectors
    II = isnan(G);                                   % Neumann indices
    if any(II), arnoldi = 0; end
    
    %----------------------------------------------------------------------
    % Generate the system matrix
    %----------------------------------------------------------------------
    n = 4*stepno;                                    % degree of polynomial term
    Np = length(pol);
    
    M = size(Z,1);
    H = zeros(n+1,n);                                % Arnoldi Hessenberg matrix
    
    %     [A, H, N] = GenSystemMatrix_Tref( ...
    %         M,H,n,Z,wc,pol,d,arnoldi,Np);
    [A, H, N, nrmlz] = GenSystemMatrix_Tref_Helm( ...
        k_wave,M,n,Z,wc,pol,d,Np);
    
    
    %----------------------------------------------------------------------
    % Associate sample points with nearest corners
    %----------------------------------------------------------------------
    Kj = zeros(M,1);
    for j = 1:M
        dd = abs(Z(j)-w);
        Kj(j) = find(dd==min(dd),1);                   % nearest corner to Zj
    end
    if rel                                            % weights to measure error
        wt = abs(Z-w(Kj))/scl;
    else
        wt = ones(M,1);
    end
    W = spdiags(sqrt(wt),0,M,M);                      % weighting for case 'rel'
    Gn = G; Gn(II) = 0;                               % set Neumann vals to 0
    
    %----------------------------------------------------------------------
    % Solve system for coefficients
    %----------------------------------------------------------------------
    c = (W*A)\(W*Gn);                                 % least-squares solution
    cc = [c(1); c(2:n+1)-1i*c(n+2:2*n+1)              % complex coeffs for f
        c(2*n+2:2*n+Np+1)-1i*c(2*n+Np+2:end)];
 
    f = @(z) fzeval_Helm(nrmlz,k,Z,wc,cc,pol,n);
    u = @(z) real(f(z));
    
    %----------------------------------------------------------------------
    % Compute error per corner
    %----------------------------------------------------------------------
    for k = 1:nw
        Kk = find(Kj==k);
        errk(k) = norm(wt(Kk).*(A(Kk,:)*c-Gn(Kk)),inf); % error near corner k
    end
    
    %----------------------------------------------------------------------
    % Global error
    %----------------------------------------------------------------------
    err = norm(wt.*(A*c-Gn),inf);                     % global error
    polmax = 100;
    
    %----------------------------------------------------------------------
    % Increment number of poles
    %----------------------------------------------------------------------
    for k = 1:nw
        if (errk(k) > q*err) & (nkv(k) < polmax)
            nkv(k) = nkv(k)+ceil(1+sqrt(nkv(k)));
        else
            nkv(k) = max(nkv(k),ceil(stepno/2));
        end
    end
    
    disp(errk)
    disp(nkv)
    
    %----------------------------------------------------------------------
    % Break if error less than prescribed tolerance
    %----------------------------------------------------------------------
    if firstFlag
        %         initErrScale =  floor(log10(err));
        %         errScale = initErrScale-2;
        %         tol = 1*10^errScale;
        %         initErrScale =  floor(log10(err));
        %         errScale = initErrScale-2;
        firstErr = err;
        tol = err*10^(-3);
        
%         if doPlot
%             figure
%             clf
%             hold on, grid on
%             plot(ww,'-k',LW,1), plot(pol,'*r',MS,6), plot(Z,'+m',MS,6),
%             set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6)
%             axis(ax), axis equal
%             xlim([-1.5, 3.5]), ylim([-1.5, 3.5])
%         end
        
        firstFlag = 0;
    end

    errvec = [errvec err]; Nvec = [Nvec; N];
    if err < .5*tol
        break
    end                        % convergence success
    
    condAvec = [condAvec cond(A)];
    coeffvec = [coeffvec norm(c)];
    
    %----------------------------------------------------------------------
    % Save the best results to date
    %----------------------------------------------------------------------
    %     if err < err0                                      % save the best so far
    %         u0 = u; f0 = f; Z0 = Z; G0 = G; A0 = A; M0 = M;
    %         N0 = N; err0 = err; pol0 = pol; wt0 = wt;
    %     end
    
    fprintf('Step: %d --- norm(g): %d \n',stepno,norm(G));
    
    figure(hFig_G)
    clf
    hold on
    plot(real(G),'b')
    plot(real(A*c),'r.')
end
warning(warn.state,'MATLAB:rankDeficientMatrix')       % back to original state
tsolve = toc;

fprintf('A: (%d,%d)\n',size(A,1),size(A,2));
fprintf('firstErr: %.16e\n',firstErr);


end   % end of main program


















% if doPlot
%     ws = 'error'; if rel, ws = 'weighted error'; end
%     hFig = figure;
%     set(hFig, 'Position', [300 100 960 680])
%     axes(PO,[.09 .65 .35 .26])
%     semilogy(sqrt(Nvec),errvec,'.-k',LW,0.7,MS,10), grid on, hold on
%     %     semilogy(sqrt(N),maxerr,'or',MS,7,LW,1), hold off
%     errmin = .01*tol;
%     axis([0 1.1*max(sqrt(Nvec)) 1e-14 100])
%     %     axis([0, 1.1*max(sqrt(Nvec)), 0.9*min(errvec), 1.1*max(errvec)])
%     set(gca,FS,fs-1), title('convergence',FS,fs,FW,NO)
%     if arnoldi == 0, title('convergence - no Arnoldi',FS,fs,FW,NO), end
%     xlabel('sqrt(DoF)',FS,fs), ylabel(ws,FS,fs)
%     set(gca,'ytick',10.^(-16:4:0))
%     ax2 = axis; x1 = ax2(1) + .05*diff(ax2(1:2));
%     s = sprintf('solve time = %6.3f secs',tsolve);
%     if ~steps, text(x1,4e-11,s,FS,fs), end
%     z = randn(1000,1)+1i*randn(1000,1); z = z/10;
%     %     tic, u(z); teval = 1e3*toc;
%     %     s = sprintf('eval time = %4.1f microsecs per pt',teval);
%     %     text(x1,4e-13,s,FS,fs)
% end
% 
% %----------------------------------------------------------------------
% % Plot solution
% %----------------------------------------------------------------------
% if doPlot
%     %     figure
%     uu = u(zz);
%     uu(~inpolygonc(zz,ww)) = nan;
%     axes(PO,[.52 .34 .47 .56]), levels = linspace(min(G),max(G),20);
%     %     contour(sx,sy,uu,levels,LW,.5), colorbar, axis equal, hold on
%     imagesc(sx,sy,uu), colorbar, axis equal, hold on, set(gca,'ydir','normal')
%     colormap(jet)
%     plot(ww,'-k',LW,1), plot(pol,'.r',MS,6)
%     set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6), axis(ax)
%     title(['dim(A) = ' int2str(M) ' x ' int2str(N) ' ', ...
%         ' #poles = ' int2str(length(pol))],FS,fs,FW,NO), hold off
% end
% 
% %----------------------------------------------------------------------
% % Plot boundary error
% %----------------------------------------------------------------------
% if doPlot
%     axes(PO,[.09 .21 .35 .28])
%     semilogy([-pi pi],maxerr*[1 1],'--b',LW,1), hold on
%     semilogy(angle(Z2-wc),wt2.*abs(u(Z2)-G2),'.r',MS,4)
%     axis([-pi pi .0001*errmin 1]), grid on
%     semilogy(angle(Z-wc),wt.*abs(u(Z)-G),'.k',MS,4), hold off
%     set(gca,'ytick',10.^(-16:4:0))
%     set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
%     set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
%     title([ws ' on boundary'],FS,fs,FW,NO)
% end
