%----------------------------------------------------------------------
% Generate the system matrix
%----------------------------------------------------------------------
function A = GenSysMatrix_MPS(k,N_2,z_samp,z_samp_weights,z_poles,z_ast)

%----------------------------------------------------------------------
% % Newman basis
W = z_samp-z_poles.';
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
W_ast_nrmlzd_pow = W_ast.^zero_to_N_2./W_ast_abs.^zero_to_N_2;

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

%----------------------------------------------------------------------
% Rescale to improve conditioning (hopefullly...)
z_samp_N = length(z_samp);

% hanRescale = 1./hanRescale_(han); hanRescale = hanRescale.^0;
% besRescale = 1./besRescale_(zero_to_N_2,k); besRescale = besRescale.^0;
% matRescale = repmat([hanRescale, hanRescale, besRescale, besRescale(2:end)], [z_samp_N 1]);
% Is the above indexing the wrong way around?
% if z_samp_weights(1) ~= 1
%     z_samp_weights = z_samp_weights.^0;
% %     z_samp_weights = z_samp_weights.^1;
%     matRescale = repmat(z_samp_weights, [1 size(A,2)]);
%     A = A .* matRescale;
% else
%     return
% end

end

function res = hanRescale_(han)
    res = sqrt(sum(abs(han).^2,1));
%     disp(max(res(:)));
end

function res = besRescale_(N_Range,k)
    r_rescale = 2;
    res = abs(besselj(N_Range, min(N_Range, k*r_rescale)));
end

