%----------------------------------------------------------------------
% Evaluate the field on a given grid
%----------------------------------------------------------------------
function r = EvalGrid_New(w,w_bndBox,k,dx,X,Y,Z,...
    a_re,a_im,b_re,b_im,...
    N_1,N_2,z_poles,z_ast)

nmax = 3e3; % default blocking size
X_N = length(X);
itcount=floor(X_N/nmax);
r=mod(X_N, nmax);
u=zeros(X_N,1);
for j=1:itcount+1
    if j<=itcount
        indrange=(j-1)*nmax+1:j*nmax;
    elseif r>0
        indrange=(j-1)*nmax:X_N;
    else
        break
    end
%     if ~isempty(p.nx)
%         ptemp=pointset(p.x(indrange),p.nx(indrange));
%     else
        X_temp = X(indrange);
%     end
    ut = EvalPoints(w,Z,X_temp, dx, w_bndBox ,nmax); % recursive!(depth 1)
    u(indrange) = ut;
end

% ptemp --- X
% o --- dx, bb, nmax


r = NaN*zeros(size(p.x));
      


% r = zeros(size(Z));
% 
% for j=1:N_1
%     w = Z-z_poles(j);
%     r = r + a_re(j)*newman_re_(k,w) + a_im(j)*newman_im_(k,w);
% end
% 
% w_ast = Z-z_ast;
% for j=0:N_2
%     r = r + b_re(j+1)*runge_re_(k,w_ast,j);
% end
% 
% for j=1:N_2
%     r = r + b_im(j)*runge_im_(k,w_ast,j);
% end


end

