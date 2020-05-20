function [A, H, N] = GenSystemMatrix_Tref( ...
    M,H,n,Z,wc,pol,d,arnoldi,Np)

                                % Arnoldi Hessenberg matrix
if arnoldi == 1
    Q = ones(M,1);                                % cols of Q have norm sqrt(M)
    for k = 1:n                                   % Arnoldi process
        v = (Z-wc).*Q(:,k);                        % (not yet implemented if
        for j = 1:k                                %  are Neumann BCs)
            H(j,k) = Q(:,j)'*v/M;
            v = v - H(j,k)*Q(:,j);
        end
        H(k+1,k) = norm(v)/sqrt(M);
        Q = [Q v/H(k+1,k)];
    end
else                                             % no-Arnoldi option
    Q = ((Z-wc)/scl).^(0:n);                      % (for educational purposes)
end
A = [real(Q) imag(Q(:,2:n+1))];                  % matrix for least-sq
% Dirichlet row entry pairs have the form Re (u+iv)*(a-ib) = [u v][a b]' = g
% Neumann row entry pairs have the form Im (U+iV)*(a-ib) = [V -U][a b]' = 0
% where U+iV = (u+iv)/(unit vector T in tangent direction)
% and now u and v correspond not to d/(Z-pol) but to -d/(Z-pol)^2
if Np > 0                                         % g linear => no poles
    A = [A real(d./(Z-pol)) imag(d./(Z-pol))];     % columns with poles
end

N = size(A,2);                                    % no. of cols = 2n+1+2Np
    
end