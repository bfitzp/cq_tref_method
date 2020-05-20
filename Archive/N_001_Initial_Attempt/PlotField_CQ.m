function PlotField_CQ(coeffs,p,c,...
    w,w_bndBox,z_samp,z_samp_weights,z_poles,z_ast,N_2,dx,...
    T,time_N,dt,k_threshold_Ind)

%----------------------------------------------------------------------
% Create grid
%----------------------------------------------------------------------
n = floor((w_bndBox(2)-w_bndBox(1))/dx); gx = w_bndBox(1) + dx*(0:n); % grids
n = floor((w_bndBox(4)-w_bndBox(3))/dx); gy = w_bndBox(3) + dx*(0:n);    
[X Y] = meshgrid(gx, gy);
Z = X + 1i*Y;
Z_ = Z(:);

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

coeffs = bsxfun(@times,coeffs,R.^(0:time_N));    
coeffs = fft(coeffs,[],2);

% --------------------------------------------------
% Solve Helmholtz equations over the given frequency
%   range.
% --------------------------------------------------
A = GenSysMatrix_MPS(1,N_2,Z_,[],z_poles,z_ast);
u = zeros(size(Z_,1),time_N+1);
di = zeros(size(Z_,1),time_N+1);
for l=freqs_Range
% for l=freqs_Range
    if l > k_threshold_Ind/2
        continue;
    else
        %----------------------------------------------------------------------
        % Evaluate on grid
        fprintf('Evaluating solution: Step %d of %d\n',l,k_threshold_Ind);
        
        k_l = 1i*p(R*omega^(-l))/dt;
        [u(:,l+1),di(:,l+1)] = EvalGrid_MPS(w,k_l,coeffs(:,l+1),Z_,...
            N_2,z_samp_weights,z_poles,z_ast);
    end
end

u(:,time_N+2-(1:floor(time_N/2)))=conj(u(:,2:floor(time_N/2)+1));

u = real(ifft(u,[],2));
u = bsxfun(@times,u,R.^(-(0:time_N)));

[r,c] = size(Z);
u_new = zeros(r,c,time_N+1);
di_new = zeros(r,c,time_N+1);
for l=1:time_N+1
    u_new(:,:,l) = reshape(u(:,l), size(Z));
    di_new(:,:,l) = reshape(di(:,l), size(Z));
end

u_tot = u_new;

%----------------------------------------------------------------------
% Plot fields
PlotField_MPS_CQ(gx,gy,real(u_tot),di_new(:,:,1),w_bndBox,time_N)

end