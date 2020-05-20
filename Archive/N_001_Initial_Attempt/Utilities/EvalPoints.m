function res = EvalPoints(...
    w,Z,X_temp, dx, w_bndBox ,nmax)
u = nan*zeros(size(X_temp));

ii = inside([], w,X_temp);

u(ii) = 0; 


% r = zeros(size(Z));

for j=1:N_1
    w = Z-z_poles(j);
    r = r + a_re(j)*newman_re_(k,w) + a_im(j)*newman_im_(k,w);
end

w_ast = Z-z_ast;
for j=0:N_2
    r = r + b_re(j+1)*runge_re_(k,w_ast,j);
end

for j=1:N_2
    r = r + b_im(j)*runge_im_(k,w_ast,j);
end



end

