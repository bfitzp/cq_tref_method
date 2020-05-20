%----------------------------------------------------------------------
% Evaluate the field on a given grid
%----------------------------------------------------------------------
function [u, di] = EvalGrid_MPS(w,k,coeffs,Z,...
    N_2,z_samp_weights,z_poles,z_ast)

%----------------------------------------------------------------------
% Split into blocks to reduce memory requirements
block_Max = 3e3; % default blocking size
Z_N = length(Z);

if Z_N > block_Max
    blocksFloor_N = floor(Z_N/block_Max);
    r = mod(Z_N, block_Max);
    di = zeros(Z_N,1);
    u = zeros(Z_N,1);
    for j=1:blocksFloor_N+1
        if j<=blocksFloor_N
            indrange=(j-1)*block_Max+1:j*block_Max;
        elseif r>0
            indrange=(j-1)*block_Max:Z_N;
        else
            break
        end
    %     if ~isempty(p.nx)
    %         ptemp=pointset(p.x(indrange),p.nx(indrange));
    %     else
            Z_temp = Z(indrange);
    %     end
        [ut,dit] = EvalGrid_MPS(w,k,coeffs,Z_temp,...
            N_2,z_samp_weights,z_poles,z_ast);
        di(indrange) = dit;
        u(indrange) = ut;
    end
    return
end

%----------------------------------------------------------------------
% Determine what grid points are inside the domain
di = nan*zeros(size(Z));
u = di;    
ii = inside([], w, Z);

di(ii) = 1;
u(ii) = 0;

%----------------------------------------------------------------------
% Compute the field for points inside the domain
A = GenSysMatrix_MPS(k,N_2,Z(ii),z_samp_weights,z_poles,z_ast);
u(ii) = A*coeffs;

end


