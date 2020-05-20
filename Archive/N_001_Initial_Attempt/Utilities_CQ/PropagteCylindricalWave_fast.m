function u_Bdy = PropagteCylindricalWave_fast(signal,p,c,z_0,z_samp,T,time_N,k_threshold_Ind)

% Geometry
diffs = bsxfun(@minus,z_samp,z_0);
% dist = sqrt(diffs(:,1).^2+diffs(:,2).^2);
dist = abs(diffs);
% U = @(k) 1i/4*besselh(0,1,k*dist);

% Sample signal
t = linspace(0,T,time_N+1);
g = signal(t);
dt = T/time_N;
R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

h = bsxfun(@times,g,R.^(0:time_N));    
h = fft(h,[],2);

d1 = size(1i/4*besselh(0,1,1*dist),1); 

u_Bdy = zeros(d1,time_N+1);
k_threshold = -1;
k_threshold_Ind = -1;
k_threshold_Flag = 0;
for l=freqs_Range
    k_l = 1i*p(R*omega^(-l))/dt;
    
    % If past the threshold flagged, set the field to zero
    if k_threshold_Flag
        u_Bdy(:,l+1) = zeros(size(h,1));
    % If not past the threshold flagged, set the field and check threshold
    else
        A = 1i/4*besselh(0,1,k_l*dist);
        u_Bdy(:,l+1) = A*h(:,l+1);
        if norm(u_Bdy(:,l+1),inf) < 1e-16
            u_Bdy(:,l+1) = zeros(size(h,1));
            k_threshold = k_l;
            k_threshold_Ind = l-1;
            k_threshold_Flag = 1;
        end
    end
end
u_Bdy(:,time_N+2-(1:floor(time_N/2)))=conj(u_Bdy(:,2:floor(time_N/2)+1));

u_Bdy = real(ifft(u_Bdy,[],2));
u_Bdy = bsxfun(@times,u_Bdy,R.^(-(0:time_N)));

end