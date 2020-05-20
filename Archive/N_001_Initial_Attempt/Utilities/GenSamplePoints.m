%----------------------------------------------------------------------
% Hardcoded for the moment...
%----------------------------------------------------------------------
function [z_samp,M] = GenSamplePoints(z_samp_N,...
    z_poles,z_polesCorner_N,w)

incFrac = (1:z_samp_N)/z_samp_N;

%----------------------------------------------------------------------
% Sample points - w_1
z_samp_1 = GenerateSamplePoints_(...
    z_poles(:,1),w(1),incFrac,...
    z_polesCorner_N,z_samp_N,1,1);

%----------------------------------------------------------------------
% Sample points - w_2
z_samp_2 = GenerateSamplePoints_(...
    z_poles(:,2),w(2),incFrac,...
    z_polesCorner_N,z_samp_N,1,-1);

%----------------------------------------------------------------------
% Sample points - w_3
z_samp_3 = GenerateSamplePoints_(...
    z_poles(:,3),w(3),incFrac,...
    z_polesCorner_N,z_samp_N,-1,-1);

%----------------------------------------------------------------------
% Sample points - w_4
z_samp_4 = GenerateSamplePoints_(...
    z_poles(:,4),w(4),incFrac,...
    z_polesCorner_N,z_samp_N,-1,1);

z_samp = [...
    z_samp_1;...
    z_samp_2;...
    z_samp_3;...
    z_samp_4;...
    ];

z_samp_N = length(z_samp);
M = z_samp_N;

end

function z_samp = GenerateSamplePoints_(...
    z_poles_vec,w,incFrac,...
    z_polesCorner_N,z_samp_N,a,b)
delta = abs(z_poles_vec - w);
inc = delta*incFrac;
z_samp_h = w*ones(z_polesCorner_N,z_samp_N) + 4*a*inc;
z_samp_v = w*ones(z_polesCorner_N,z_samp_N) + 4*b*1i*inc;

z_samp = [z_samp_h(:);z_samp_v(:)];

threshold = 0.99;
z_samp = z_samp(find(abs(z_samp-w) < threshold));

end