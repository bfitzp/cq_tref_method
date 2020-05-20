function [u, maxerr, f, Z, Zplot, A] = N_0002_Helmholtz_Eq_(k_wave, P, varargin)

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
wc = mean(wr+1i*wi);                             % for scale- and transl-invariance
scl = max([diff(wr),diff(wi)]);

%----------------------------------------------------------------------
% Modifications
%----------------------------------------------------------------------
z_0 = 0.5 - 0.0i;
for k = 1:nw
    g{k} = @(z) real(1i/4*besselh(0,k_wave*abs(z-z_0)));       % default
end
tol = 1e-10;
% wc = 0.8-0.5i;
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
        dk = exp(2*sk);
        dk = scl*dk;                                  % stronger clustering near corner
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
    
    %----------------------------------------------------------------------
    % Generate sample points and boundary data
    %----------------------------------------------------------------------
    for k = 1:nw
        tt{k} = sort(tt{k}(:));
        tk = tt{k}; pk = pt{k};                       % abbrevations
        Z = [Z; pk(tk)];                              % sample pts on side k
        G = [G; g{k}(pk(tk))];                        % boundary data at these pts
        h = 1e-4;                                     % 4-pt trapezoidal rule
        T = [T; (pk(tk+h)-1i*pk(tk+1i*h) ...
            - pk(tk-h)+1i*pk(tk-1i*h))/(4*h);];    % unnormalized tangent vectors
    end
    T = T./abs(T);                                   % normalize tangent vectors
    II = isnan(G);                                   % Neumann indices
    if any(II), arnoldi = 0; end
    
    %----------------------------------------------------------------------
    % Generate the system matrix
    %----------------------------------------------------------------------
    n = 4*stepno;                                    % degree of polynomial term
    Np = length(pol);
    
    M = size(Z,1);
    H = zeros(n+1,n);                                % Arnoldi Hessenberg matrix
    
    if stepno == 50
        aaa = 1;
    end
    %     [A, H, N] = GenSystemMatrix_Tref( ...
    %         M,H,n,Z,wc,pol,d,arnoldi,Np);
    [A, H, N, nrmlz] = GenSystemMatrix_Tref_Helm( ...
        k_wave,M,H,n,Z,wc,pol,d,arnoldi,Np);
    
    
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
    
    % %     %----------------------------------------------------------------------
    % %     % Declare evaluation function and solution function
    % %     %----------------------------------------------------------------------
    %     f = @(z) reshape(fzeval(z(:),wc,...               % vector and matrix inputs
    %         cc,H,pol,d,arnoldi,scl,n),size(z));    % to u and f both allowed
    %     u = @(z) real(f(z));
    
    f = @(z) reshape(fzeval_Helm(nrmlz,k_wave,z(:),wc,...               % vector and matrix inputs
        c,H,pol,d,arnoldi,scl,n),size(z));    % to u and f both allowed
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
    %     polmax = 100;
    polmax = 100;
    
    %----------------------------------------------------------------------
    % Increment number of poles
    %----------------------------------------------------------------------
    for k = 1:nw
        if (errk(k) > q*err) & (nkv(k) < polmax)
            %             if k==4
            nkv(k) = nkv(k)+ceil(1+sqrt(nkv(k)));      % increase no. poles
            %                 nkv(k) = nkv(k)+ceil(1+(nkv(k).^(1/4)));      % increase no. poles
            %             else
            %                 nkv(k) = max(nkv(k),ceil(stepno/2));
            %             end
        else
            nkv(k) = max(nkv(k),ceil(stepno/2));
            %nkv(k) = min(polmax,nkv(k)+1);
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
        tol = err*10^(-8);
        figure
        clf
        hold on, grid on
        plot(ww,'-k',LW,1), plot(pol,'*r',MS,6), plot(Z,'+m',MS,6),
        set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6)
        axis(ax), axis equal
        xlim([-1.5, 3.5]), ylim([-1.5, 3.5])
        %        title(['dim(A) = ' int2str(M) ' x ' int2str(N) ' ', ...
        %            ' #poles = ' int2str(length(pol))],FS,fs,FW,NO), hold off
        firstFlag = 0;
    end
    
    
    errvec = [errvec err]; Nvec = [Nvec; N];
    if err < .5*tol, break, end                        % convergence success
    
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
end
warning(warn.state,'MATLAB:rankDeficientMatrix')       % back to original state
tsolve = toc;  % =========== end of main loop =========================================

fprintf('A: (%d,%d)\n',size(A,1),size(A,2));
% geomFields_FigPos = [900 100 760 900];
% z_poles = pol;
% N_2 = n;
% z_ast = wc;
% coeffs = c;
% dx = 0.01;
%
% figure
% Plot_Geom_and_Fields(...
%     geomFields_FigPos,dx,coeffs,w,w_bndBox,...
%     k_wave,z_poles.',Z,[],z_ast,N_2);

%----------------------------------------------------------------------
% A posteriori error check using finer mesh
%----------------------------------------------------------------------
[err2, maxerr, Z2, G2, wt2] = APosterioriErrorCheck(nw,tt,pt,g,w,scl,u,rel,err);

%----------------------------------------------------------------------
% Plot convergence
%----------------------------------------------------------------------
% figure
% % semilogy(sqrt(Nvec),condAvec), grid on, hold on
% semilogy(sqrt(Nvec),coeffvec), grid on, hold on
% xlim([0 1.1*max(sqrt(Nvec))])
%
if plots
    ws = 'error'; if rel, ws = 'weighted error'; end
    hFig = figure;
    set(hFig, 'Position', [100 50 860 680])
    %     if steps
    %         figure
    %     else
    %         clf
    %     end
    %     shg
    axes(PO,[.09 .65 .35 .26])
    semilogy(sqrt(Nvec),errvec,'.-k',LW,0.7,MS,10), grid on, hold on
    %     semilogy(sqrt(N),maxerr,'or',MS,7,LW,1), hold off
    errmin = .01*tol;
    axis([0 1.1*max(sqrt(Nvec)) 1e-14 100])
    %     axis([0, 1.1*max(sqrt(Nvec)), 0.9*min(errvec), 1.1*max(errvec)])
    set(gca,FS,fs-1), title('convergence',FS,fs,FW,NO)
    if arnoldi == 0, title('convergence - no Arnoldi',FS,fs,FW,NO), end
    xlabel('sqrt(DoF)',FS,fs), ylabel(ws,FS,fs)
    set(gca,'ytick',10.^(-16:4:0))
    ax2 = axis; x1 = ax2(1) + .05*diff(ax2(1:2));
    s = sprintf('solve time = %6.3f secs',tsolve);
    if ~steps, text(x1,4e-11,s,FS,fs), end
    z = randn(1000,1)+1i*randn(1000,1); z = z/10;
    %     tic, u(z); teval = 1e3*toc;
    %     s = sprintf('eval time = %4.1f microsecs per pt',teval);
    %     text(x1,4e-13,s,FS,fs)
end

%----------------------------------------------------------------------
% Plot solution
%----------------------------------------------------------------------
if plots
    %     figure
    uu = u(zz); uu(~inpolygonc(zz,ww)) = nan;
    axes(PO,[.52 .34 .47 .56]), levels = linspace(min(G),max(G),20);
    %     contour(sx,sy,uu,levels,LW,.5), colorbar, axis equal, hold on
    imagesc(sx,sy,uu), colorbar, axis equal, hold on, set(gca,'ydir','normal')
    colormap(jet)
    plot(ww,'-k',LW,1), plot(pol,'.r',MS,6)
    set(gca,FS,fs-1), plot(real(wc),imag(wc),'.k',MS,6), axis(ax)
    title(['dim(A) = ' int2str(M) ' x ' int2str(N) ' ', ...
        ' #poles = ' int2str(length(pol))],FS,fs,FW,NO), hold off
end

%----------------------------------------------------------------------
% Plot boundary error
%----------------------------------------------------------------------
if plots
    axes(PO,[.09 .21 .35 .28])
    semilogy([-pi pi],maxerr*[1 1],'--b',LW,1), hold on
    semilogy(angle(Z2-wc),wt2.*abs(u(Z2)-G2),'.r',MS,4)
    axis([-pi pi .0001*errmin 1]), grid on
    semilogy(angle(Z-wc),wt.*abs(u(Z)-G),'.k',MS,4), hold off
    set(gca,'ytick',10.^(-16:4:0))
    set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
    set(gca,FS,fs-1), xlabel('angle on boundary wrt wc',FS,fs)
    title([ws ' on boundary'],FS,fs,FW,NO)
end

end   % end of main program





% function fZ = fzeval(Z,wc,cc,H,pol,d,arnoldi,scl,n)
% ZZ = [wc; Z];
% if arnoldi
%     Q = ones(size(ZZ));
%     for k = 1:size(H,2)
%         v = (ZZ-wc).*Q(:,k);
%         for j = 1:k
%             v = v - H(j,k)*Q(:,j);
%         end
%         Q = [Q v/H(k+1,k)];
%     end
% else
%     Q = ((ZZ-wc)/scl).^(0:n);
% end
% if length(pol) > 0
%     fZZ = [Q d./(ZZ-pol)]*cc;
% else
%     fZZ = Q*cc;
% end
% fZ = fZZ(2:end) - 1i*imag(fZZ(1));
% end














%LAPLACE  Lightning Laplace solver.
%         U = LAPLACE(P,G) solves the Laplace equation with Dirichlet or
%         homogeneous Neumann boundary data on the simply-connected region
%         Omega bounded by P, which may be a polygon or circular polygon,
%         with computation also of the harmonic conjugate.
%         v6, (c) Lloyd N. Trefethen, U. of Oxford, March 2020.
%         For info see https://people.maths.ox.ac.uk/trefethen/lightning.html.
%         Please send corrections or suggestions to trefethen@maths.ox.ac.uk.
%
%  Inputs:
%      P = vector of corners as complex numbers z = x+iy in counterclockwise
%              order to specify a polygon
%          or cell array of corners v and pairs [v r] to specify a circular
%              polygon: r = radius of curvature of arc from this v to the next
%          or 'pent'[agon], 'snow'[flake], 'iso'[spectral], 'L', or 'circleL'
%          or integer >= 3, the number of corners of a random polygon
%          or integer <= -3, -1 x no. of corners of a random circular polygon
%
%      g = function handle for real-valued Dirichlet boundary data
%          or cell array of function handles for sides P1-P2, P2-P3,...
%          or vector of constant values for these sides
%          (default @(z) real(z).^2)).  If g = nan on any side, the
%          BC imposed there is homogeneous Neumann, i.e., du/dn = 0.
%
%  Further optional inputs:
%    'tol', tol = tolerance for maximal absolute error (default 1e-6)
%    'noplots' to suppress plotting
%    'steps' for step-by-step plots of errors on boundary and poles
%    'rel' to weight error by scaled dist. to corner (automatic if g discont.)
%    'slow' to turn off adaptive mode for cleaner root-exp convergence curves
%    'noarnoldi' to turn off Arnoldi stabilization of polynomial term
%    'aaa' for AAA compression of result.  Requires Chebfun in path.
%
%  Outputs for [U,MAXERR,F,Z,ZPLOT,A] = laplace(P,G)
%       u = function handle for solution u(z) of lap(u) = 0, u = g on boundary
%  maxerr = upper bound on maximal error, even near singularities at corners
%             (or error weighted by distances to corners, if 'rel' is specified)
%       f = function handle for analytic function f = u + iv
%       Z = sample points on boundary
%       Zplot = set of points for plotting boundary
%       A = final rectangular matrix used for least-squares problem
%
% Reference: A. Gopal and L. N. Trefethen, Solving Laplace problems with corner
% singularities via rational functions, SIAM J Numer Anal 57 (2019), 2074-94.
%
% Examples:
%
%   laplace([0 1 1+1i 1i],[0 0 0 1]);       % square with piecewise const bc
%   laplace('iso',1:8);                     % isospectral octagon
%   laplace('iso',[0 nan 1 1 1 nan 0 0]);   % same but with some Neumann BCs
%   laplace(8,@(z) exp(real(z)));           % random octagon with BC exp(x)
%   laplace(-8,@(z) exp(real(z)));          % same but circular arc boundaries
%   laplace({[1 .5] 1+1i -1+1i -1});        % bullet
%   laplace('circleL',[0 1 0 0 0]);         % circular L-shape
%   [u,maxerr] = laplace('L','tol',1e-10);  % see NA Digest, Nov. 2018
%       u(.99+.99i), maxerr                 %    exact soln. 1.0267919261073...
