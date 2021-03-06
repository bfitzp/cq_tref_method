function [coeffs, coeffs_fft, nrmlz_vec, n, u_Bdy_coeffs, testFreqResults] = SolveCQ( ...
    P,u_Bdy,u_Bdy_fft,p,c, ...
    Z,pol,d, ...
    T,time_N,dt,k_threshold_Ind, ...
    tt,pt,g, ...
    ww,wc,nw, ...
    z_0, ...
    Z_hr, ...
    u_in_hr,u_in_fft_hr)

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

% b = bsxfun(@times,u_Bdy,R.^(0:time_N));
% b = fft(b,[],2);
b = u_Bdy_fft;
b_hr = u_in_fft_hr;

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
coeffs = zeros(size(A,2),time_N+1);
nrmlz_vec = zeros(size(A,2),time_N+1);

testFreqResults = {};
% parfor l=freqs_Range
for l=freqs_Range
    k_l = 1i*p(R*omega^(-l))/(c*dt);
    % Adaptively solve the first Helmholtz equation
    if l==0
        [c_adapt, nrmlz_adapt, Z_adapt, pol_adapt, ...
            scl_adapt, n_adapt, d_adapt] = SolveCQ_Adapt(...
            k_l, P, z_0, Z_hr, b_hr(:,l+1));
        
        testFreqResults = {c_adapt, nrmlz_adapt, Z_adapt, pol_adapt, ...
            scl_adapt, n_adapt, d_adapt};
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

        coeffs(:,l+1) = A\b(:,l+1);
%         coeffs(:,l+1) = lsqminnorm(A,b(:,l+1));
%         coeffs = pinv(A)*b(:,l+1);    
    end
end

coeffs_fft = coeffs;

coeffs(:,time_N+2-(1:floor(time_N/2))) = conj(coeffs(:,2:floor(time_N/2)+1));
% u_Bdy_coeffs(:,time_N+2-(1:floor(time_N/2))) = conj(u_Bdy_coeffs(:,2:floor(time_N/2)+1));

coeffs = ifft(coeffs,[],2);
coeffs = bsxfun(@times,coeffs,R.^(-(0:time_N)));

end