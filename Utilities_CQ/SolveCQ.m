function testFreqResults = SolveCQ( ...
    P, ...
    u_in_bdy_fft_hr, ...
    p,c, ...
    time_N,dt, ...
    z_0, ...
    Z_hr, ...
    err_tol)

%----------------------------------------------------------------------
% Parameters needed for calculating Laplace parameter
%----------------------------------------------------------------------
R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

u_in_bdy_hr_ = ifft(u_in_bdy_fft_hr,[],2);
err_tol_ast = err_tol* norm(u_in_bdy_hr_(:),inf);

% --------------------------------------------------
% Solve Helmholtz equations over the given frequency
%   range.
% --------------------------------------------------
% maxstepno = 30;
% M = size(Z,1);
% n = 4*maxstepno;                                    % degree of polynomial term
% Np = length(pol);

%----------------------------------------------------------------------
% Adaptively solve Helmholtz equations
%----------------------------------------------------------------------
testFreqResults = cell(1);
l = 0;
tolAchieved = 0;
while ~tolAchieved
    fprintf('Solving. l = %d\n',l);
    k_l = 1i*p(R*omega^(-l))/(c*dt);
    
    [finishedFlag, c_adapt, nrmlz_adapt, pol_adapt, ...
        n_adapt] = SolveHelmholtzAdaptive(...
        k_l, P, Z_hr, -u_in_bdy_fft_hr(:,l+1), err_tol_ast);
        
    if finishedFlag == 0
        testFreqResults{l+1} = {c_adapt, nrmlz_adapt, pol_adapt, ...
            n_adapt};
    else
        tolAchieved = 1;
    end
    l = l+1;
end

end