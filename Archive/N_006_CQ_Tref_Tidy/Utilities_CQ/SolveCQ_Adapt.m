function [c_adapt, nrmlz_adapt, Z_adapt, pol_adapt, ...
    n_adapt, d_adapt] = SolveCQ_Adapt(...
    k, P, z_0, Z_, b_)

%----------------------------------------------------------------------
% Solve this Helmholtz equation adaptively with Trefethen's method
%----------------------------------------------------------------------
% For this wavenumber, this will give us a new set of each of the following:
%   1) coefficients
%   2) normalization vector
%   3) Sample points
%   4) Poles
%----------------------------------------------------------------------
[c_adapt, nrmlz_adapt, Z_adapt, pol_adapt, ...
    n_adapt, d_adapt] = N_0003_Helmholtz_Eq_CQ_Export_(...
    k, P, z_0, Z_, b_);

end

