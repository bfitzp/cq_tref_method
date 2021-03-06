function [Z,pol,d,tt,pt,g, ...
    ww,wr,wi,wc,nw,scl] = GenSamplePointsAndBdyData(...
    z_0, P, poles_N, varargin)

%% Set up the problem
[g, w, ww, pt, dw, tol, steps, plots, ...        % parse inputs
    slow, rel, arnoldi, aaaflag] = ...
    parseinputs(P,varargin{:});
Zplot = ww;
nw = length(w);                                  % number of corners
wr = sort(real(ww)); wr = wr([1 end]);
wi = sort(imag(ww)); wi = wi([1 end]);
wc = mean(wr+1i*wi);                             % for scale- and transl-invariance
scl = max([diff(wr),diff(wi)]);

% %----------------------------------------------------------------------
% % Incident field
% %----------------------------------------------------------------------
for k = 1:nw
    g{k} = @(k_wave,z) real(1i/4*besselh(0,k_wave*abs(z-z_0)));       % default
end

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

%----------------------------------------------------------------------
% Main loop: increase number of poles until convergence
%----------------------------------------------------------------------
Nvec = []; errvec = []; condAvec = []; coeffvec = [];
tic
errk = ones(nw,1);                               % max error near each corner
nkv = zeros(nw,1);                               % no. of poles at each corner
minstepno = 1;
maxstepno = 35; err0 = Inf;

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
    nk = poles_N;                                  % no. of poles at this corner
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
        linspace(0,dw(k),nk)];            % additional pts along side
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
%     G = [G; g{k}(pk(tk))];                        % boundary data at these pts
%     h = 1e-4;                                     % 4-pt trapezoidal rule
%     T = [T; (pk(tk+h)-1i*pk(tk+1i*h) ...
%         - pk(tk-h)+1i*pk(tk-1i*h))/(4*h);];    % unnormalized tangent vectors
end
% T = T./abs(T);                                   % normalize tangent vectors
% II = isnan(G);                                   % Neumann indices
% if any(II), arnoldi = 0; end
  


end