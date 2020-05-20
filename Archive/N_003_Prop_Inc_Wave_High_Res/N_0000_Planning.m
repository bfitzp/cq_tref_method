%----------------------------------------------------------------------
% Generalizing Trefethen's Method to the time domain
%----------------------------------------------------------------------
% Propagate incident wave to the boundary
%----------------------------------------------------------------------
%   Upsample the time resolution - ON HOLD
%   Use a very fine spatial grid so we can obtain accurate values of the
%       incident field at various frequencies through interpolation
%   We set end time T and number of timesteps N
%       This gives the time grid and frequency range we need to solve over
%   We must choose a maximal set of sample points
%       1) Uniform points - 
%           Find max oscillation out of incident field
%               Eg. k_re = 4
%               Use Nyquist sampling -> Samp frequency = 8
%                   See 3D Seismic imaing by Biondi, p. 104
%       2) Near corner samp points --- Use max number of poles
%       Choose maximum 
%----------------------------------------------------------------------
% Next task:
%   1. Set T and N
%       Upsample timesteps for finer resolution
%   2. Get frequency range
%   3. Get maximum frequency
%   4. Calculate max samp frequency
%   5. Compute uniform samp points --- NEED SYNC
%       Plot geometry
%           Fix for squares
%           Fix for Lshape
%   6. Use max poles to compute near corner samp points
%   7. Merge the list of points
%   8. Propagte the incident wave to the samp points
%   9. Back to time domain