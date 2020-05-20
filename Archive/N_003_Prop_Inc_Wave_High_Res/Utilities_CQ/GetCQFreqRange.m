function [k_Range, k_Nyq] = GetCQFreqRange(time_N,dt,p)
R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

k_Range = zeros(floor((time_N+1)/2),1);
for j = 1:floor((time_N+1)/2)
    k_j = 1i*p(R*omega^(-j))/dt;
    k_Range(j) = k_j;
end

k_Nyq = ceil(2 * max(abs(real(k_Range))));

end

% [m,i] = max(abs(real(k_Range)));
% k_ = k_Range(i);
% k_re_ = real(k_Range(i));
% 
% 
% 
% x = linspace(-pi,pi,200);
% % f_x_1 = real(exp(1i*k_re_*x));
% f_x_1 = real(exp(1i*(k_re_+0i)*x));
% 
% samp_pts_0 = x(1) + (0:k_Nyq)/k_Nyq * 2*pi;
% 
% figure
% hold on, grid on
% plot(x,f_x_1,'b')
% plot(samp_pts_0,0*samp_pts_0,'r*')
% xlim([min(x),max(x)])
