function res = newman_im_(k,w)
res = besselh(1,k*abs(w)).*imag(w./abs(w));
% res = besselh(1,k*abs(w)).*(real(w./abs(w)) + imag(w./abs(w)));
end

