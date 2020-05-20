function res = runge_re_(k,w,j)
res = besselj(j,k*abs(w)).*real(w.^j./abs(w.^j));
end

