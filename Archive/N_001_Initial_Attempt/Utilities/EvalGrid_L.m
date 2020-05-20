%----------------------------------------------------------------------
% Evaluate the field on a given grid
%----------------------------------------------------------------------
function r = EvalGrid_L(k,Z,...
    a_re,a_im,b_re,b_im,...
    N_1,N_2,z_poles,z_ast)

r = zeros(size(Z));

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


for i=1:size(Z,1)
    for j=1:size(Z,2)
        if real(Z(i,j)) < 0 && imag(Z(i,j)) > 0
            r(i,j) = nan;
        end
    end
end

end

