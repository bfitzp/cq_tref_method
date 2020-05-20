function [coeffs, coeffs_fft, AScale] = SolveCQ(u_Bdy,p,c,...
    z_samp,z_samp_weights,z_poles,z_ast,N_2,...
    T,time_N,dt,k_threshold_Ind)

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

b = bsxfun(@times,u_Bdy,R.^(0:time_N));    
b = fft(b,[],2);

% --------------------------------------------------
% Solve Helmholtz equations over the given frequency
%   range.
% --------------------------------------------------
A = GenSysMatrix_MPS(1,N_2,z_samp,z_samp_weights,z_poles,z_ast);
coeffs = zeros(size(A,2),time_N+1);
AScale = zeros(1,time_N+1);
parfor l=freqs_Range
    if l > k_threshold_Ind/1
        continue;
    else
        fprintf('Solving: Step %d of %d\n',l,k_threshold_Ind);
        
        k_l = 1i*p(R*omega^(-l))/dt;  
        A = GenSysMatrix_MPS(k_l,N_2,z_samp,z_samp_weights,z_poles,z_ast);
        
        AScale_ = 1/norm(A);
        AScale_ = 1;
        A = AScale_*A;
        AScale(:,l+1) = AScale_;
        
        coeffs(:,l+1) = A\b(:,l+1);
    end
end

coeffs(:,time_N+2-(1:floor(time_N/2))) = conj(coeffs(:,2:floor(time_N/2)+1));

coeffs_fft = coeffs;

coeffs = real(ifft(coeffs,[],2));
coeffs = bsxfun(@times,coeffs,R.^(-(0:time_N)));

end