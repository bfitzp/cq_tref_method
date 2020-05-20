function fZ = fzeval_Helm(nrmlz,k,Z,wc,cc,H,pol,d,arnoldi,scl,n)

z_samp = Z;
z_poles = pol;
N_2 = n;
z_ast = wc;

%----------------------------------------------------------------------
% % Newman basis
W = z_samp-z_poles;
W_abs = abs(W);
W_nrmlzd_pow = W./W_abs;

r = k*W_abs;
han = besselh(1, r);

A_1 = han.*real(W_nrmlzd_pow);
A_2 = han.*imag(W_nrmlzd_pow);

%----------------------------------------------------------------------
% Runge basis
zero_to_N_2 = 0:N_2;

W_ast = z_samp-z_ast;
W_ast_abs = abs(W_ast);
% W_ast_nrmlzd_pow = W_ast.^zero_to_N_2./W_ast_abs.^zero_to_N_2;
W_ast_nrmlzd_pow = (W_ast./W_ast_abs).^zero_to_N_2;

ii = find(W_ast==0);
W_ast_nrmlzd_pow(ii,:) = 1;

r = k*W_ast_abs;
bes = besselj(repmat(zero_to_N_2, [numel(r) 1]), repmat(r(:), [1, N_2+1]));

A_3 = bes.*real(W_ast_nrmlzd_pow);
A_4 = bes(:,2:end).*imag(W_ast_nrmlzd_pow(:,2:end));

%----------------------------------------------------------------------
% System matrix
A = [...
    A_1,...
    A_2,...
    A_3,....
    A_4...
    ];

A = bsxfun(@rdivide,A,nrmlz);

fZ = A*cc;

end