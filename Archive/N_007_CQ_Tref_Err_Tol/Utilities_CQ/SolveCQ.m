function [coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs, testFreqResults] = SolveCQ( ...
    P,u_in_bdy_fft,p,c, ...
    Z,pol,d, ...
    T,time_N,dt,k_threshold_Ind, ...
    tt,pt,g, ...
    ww,wc,nw, ...
    z_0, ...
    Z_hr, ...
    u_in_bdy_fft_hr, ...
    err_tol)

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

b = u_in_bdy_fft;
b_hr = u_in_bdy_fft_hr;

u_in_bdy_ = ifft(b,[],2);
u_in_bdy_hr_ = ifft(b_hr,[],2);
err_tol_ast = err_tol* norm(u_in_bdy_(:),inf);
% --------------------------------------------------
% Solve Helmholtz equations over the given frequency
%   range.
% --------------------------------------------------
maxstepno = 30;
M = size(Z,1);
n = 4*maxstepno;                                    % degree of polynomial term
Np = length(pol);

[A, ~, ~, ~] = GenSystemMatrix_Tref_Helm( ...
    1,M,n,Z,wc,pol,d,Np);
u_Bdy_coeffs = zeros(length(b),time_N+1);
coeffs_fft = zeros(size(A,2),time_N+1);
nrmlz_vec = zeros(size(A,2),time_N+1);

testFreqResults = cell(2,1);
doAdaptive = 1;
% parfor l=freqs_Range
for l=freqs_Range
    k_l = 1i*p(R*omega^(-l))/(c*dt);
    % Adaptively solve the first Helmholtz equation
    if doAdaptive
        [finishedFlag, c_adapt, nrmlz_adapt, Z_adapt, pol_adapt, ...
            n_adapt, d_adapt] = N_0003_Helmholtz_Eq_CQ_Export_(...
            k_l, P, z_0, Z_hr, b_hr(:,l+1), err_tol_ast);
        
        if finishedFlag == 0
            testFreqResults{l+1,1} = {c_adapt, nrmlz_adapt, Z_adapt, pol_adapt, ...
                n_adapt, d_adapt};
        else
            doAdaptive = 0;
        end
    end
    if l > k_threshold_Ind
        continue;
    else
        fprintf('Solving: Step %d of %d\n',l,k_threshold_Ind);

        %----------------------------------------------------------------------
        % Generate the system matrix
        %----------------------------------------------------------------------
        [A, ~, ~, nrmlz] = GenSystemMatrix_Tref_Helm( ...
            k_l,M,n,Z,wc,pol,d,Np);
        nrmlz_vec(:,l+1) = nrmlz;

        coeffs_fft(:,l+1) = A\b(:,l+1);
%         coeffs(:,l+1) = lsqminnorm(A,b(:,l+1));
%         coeffs = pinv(A)*b(:,l+1);    
    end
end

% coeffs(:,time_N+2-(1:floor(time_N/2))) = conj(coeffs(:,2:floor(time_N/2)+1));
% % u_Bdy_coeffs(:,time_N+2-(1:floor(time_N/2))) = conj(u_Bdy_coeffs(:,2:floor(time_N/2)+1));
% 
% coeffs = ifft(coeffs,[],2);
% coeffs = bsxfun(@times,coeffs,R.^(-(0:time_N)));

end