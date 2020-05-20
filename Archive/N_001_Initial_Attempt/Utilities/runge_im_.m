function res = runge_im_(k,w,j)
res = besselj(j,k*abs(w)).*imag(w.^j./abs(w.^j));
end

