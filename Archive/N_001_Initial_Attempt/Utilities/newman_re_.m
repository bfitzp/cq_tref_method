function res = newman_re_(k,w)
res = besselh(1,k*abs(w)).*real(w./abs(w));
% res = besselh(0,k*abs(w));
end

