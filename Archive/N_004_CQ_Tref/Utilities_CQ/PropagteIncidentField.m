function [u_Bdy, u_Bdy_fft, k_threshold_Ind] =...
    PropagteIncidentField(...
    signal,p,c,z_0,Z,T,time_N)

% Spatial information
diffs = bsxfun(@minus,Z,z_0);
dist = abs(diffs);

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

u_Bdy_norm_first = -1;
k_threshold_Ind = -1;
for l=freqs_Range

    k_l = 1i*p(R*omega^(-l))/(c*dt);
    A = 1i/4*besselh(0,1,k_l*dist);
    u_Bdy(:,l+1) = A*h(:,l+1);

%     fprintf('l: %d   ---   k_l: %.4e i%.4e\n', l, real(k_l), imag(k_l));
    fprintf('l: %d   ---   k_l:   %.2e i%.2e   ---   norm(h): %.4e   ---   norm(A): %.4e   ---   norm(u_in): %.4e\n', ...
            l, real(k_l), imag(k_l), norm(h(:,l+1)),norm(A),norm(u_Bdy(:,l+1)));

    if u_Bdy_norm_first == -1
        u_Bdy_norm_first = norm(u_Bdy(:,l+1));
%     end
    elseif norm(u_Bdy(:,l+1)) < u_Bdy_norm_first*1e-4
% %     if l == 10
%         u_Bdy(:,l+1) = u_Bdy(:,l+1).*0;
        k_threshold_Ind = l;
        break;
    end
            
end
u_Bdy(:,time_N+2-(1:floor(time_N/2))) = conj(u_Bdy(:,2:floor(time_N/2)+1));
u_Bdy_fft = u_Bdy;

u_Bdy = ifft(u_Bdy,[],2);

u_Bdy = bsxfun(@times,u_Bdy,R.^(-(0:time_N)));

end