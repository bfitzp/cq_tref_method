function [k_threshold,k_threshold_Ind,k_re_max] = ...
    ComputeFreqThreshold(threshold,time_N,dt,p,z_samp,incFieldOpt,plotFlag)

R = eps^(0.5/(time_N+1));
omega = exp(2*pi*1i/(time_N+1));

freqs_N = floor((time_N+1)/2);
freqs_Range = 0:freqs_N;

bdyNorms = zeros(freqs_N+1,1);
k_threshold = -1;
k_threshold_Ind = -1;
k_threshold_Flag = 0;
k_re_max = -1;
for l = freqs_Range
    k_l = 1i*p(R*omega.^(-l))/dt;
    u_in_ = SelectIncidentField(k_l,incFieldOpt);
    b = u_in_(z_samp);
    bdyNorms(l+1) = norm(b);
    if ~k_threshold_Flag
        if bdyNorms(l+1) < threshold
            k_threshold = 1i*p(R*omega.^(-(l-1)))/dt;
            k_threshold_Ind = l-1;
            k_threshold_Flag = 1;
        end
    end
    if abs(real(k_l)) > k_re_max
        k_re_max = abs(real(k_l));
    end
end

if plotFlag
    figure
    hold on, grid on
    plot(freqs_Range,bdyNorms)
    plot([0,freqs_N],[1e-16,1e-16])
    xlim([0,freqs_N])
    set(gca,'yscale','log')
end

end

