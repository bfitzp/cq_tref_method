function PlotScatteredFieldOnBdy_CQ(...
    coeffs,coeffs_fft,w,z_samp,z_poles,z_ast,...
    cornerInds,N_2,...
    p,time_N,dt,k_threshold_Ind,AScale)

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

% coeffs = bsxfun(@times,coeffs,R.^(0:time_N));    
% coeffs = fft(coeffs,[],2);

% --------------------------------------------------
% Evaluate on grid
% --------------------------------------------------
A = GenSysMatrix_MPS(1,N_2,z_samp,[],z_poles,z_ast);
u = zeros(size(z_samp,1),time_N+1);
parfor l=freqs_Range
% for l=freqs_Range
    if l > k_threshold_Ind/1
        continue;
    else
        fprintf('Evaluating solution: Step %d of %d\n',l,k_threshold_Ind);
        
        k_l = 1i*p(R*omega^(-l))/dt;
        u(:,l+1) = EvalGrid_MPS(w,k_l,AScale(l+1)*coeffs_fft(:,l+1),z_samp,...
            N_2,[],z_poles,z_ast);
    end
end

u(:,time_N+2-(1:floor(time_N/2))) = conj(u(:,2:floor(time_N/2)+1));

u = real(ifft(u,[],2));
u = bsxfun(@times,u,R.^(-(0:time_N)));

Plot_BdyField_CQ(z_samp,cornerInds,real(u),1,time_N,dt);

end