%----------------------------------------------------------------------
% Hardcoded for the moment...
%----------------------------------------------------------------------
function [z_poles,N_1,w_vec]...
    = GenPoles(...
    z_polesCorner_N,z_polesOffset,sigma,...
    w_1,w_2,w_3,w_4)

w_vec = [w_1;w_2;w_3;w_4];
corners_N = length(w_vec);

j = (0:z_polesCorner_N-1).';

% scale = 1;
offset_1 = z_polesOffset*(w_1);
offset_2 = z_polesOffset*(w_2);
offset_3 = z_polesOffset*(w_3);
offset_4 = z_polesOffset*(w_4);
% z_poles_1 = offset_1 - exp(-sigma*j/sqrt(z_polesCorner_N))*exp(1i*pi/4); z_poles_1 = scale*z_poles_1 + w_1;
% z_poles_2 = offset_2 - exp(-sigma*j/sqrt(z_polesCorner_N))*exp(-1i*pi/4); z_poles_2 = scale*z_poles_2 + w_2;
% z_poles_3 = offset_3 + exp(-sigma*j/sqrt(z_polesCorner_N))*exp(1i*pi/4); z_poles_3 = scale*z_poles_3 + w_3;
% z_poles_4 = offset_4 + exp(-sigma*j/sqrt(z_polesCorner_N))*exp(-1i*pi/4); z_poles_4 = scale*z_poles_4 + w_4;
% z_poles = [z_poles_1,z_poles_2,z_poles_3,z_poles_4];

a_2 = w_2-w_1;
a_1 = w_4-w_1;
c = norm(a_2)*a_1 + norm(a_1)*a_2;
c_nrmlzd = c/norm(c);

z_poles_1 = w_1 - exp(-sigma*j/sqrt(z_polesCorner_N))*c_nrmlzd;

a_2 = w_3-w_2;
a_1 = w_1-w_2;
c = norm(a_2)*a_1 + norm(a_1)*a_2;
c_nrmlzd = c/norm(c);

z_poles_2 = w_2 - exp(-sigma*j/sqrt(z_polesCorner_N))*c_nrmlzd;

a_2 = w_4-w_3;
a_1 = w_2-w_3;
c = norm(a_2)*a_1 + norm(a_1)*a_2;
c_nrmlzd = c/norm(c);

z_poles_3 = w_3 - exp(-sigma*j/sqrt(z_polesCorner_N))*c_nrmlzd;

a_2 = w_1-w_4;
a_1 = w_3-w_4;
c = norm(a_2)*a_1 + norm(a_1)*a_2;
c_nrmlzd = c/norm(c);

z_poles_4 = w_4 - exp(-sigma*j/sqrt(z_polesCorner_N))*c_nrmlzd;

z_poles_1 = z_poles_1 + offset_1;
z_poles_2 = z_poles_2 + offset_2;
z_poles_3 = z_poles_3 + offset_3;
z_poles_4 = z_poles_4 + offset_4;

z_poles = [z_poles_1,z_poles_2,z_poles_3,z_poles_4];

N_1 = corners_N*z_polesCorner_N;

% figure
% hold on, grid on
% scatter(real(z_poles_1),imag(z_poles_1));
% scatter(real(z_poles_2),imag(z_poles_2));
% scatter(real(z_poles_3),imag(z_poles_3));
% scatter(real(z_poles_4),imag(z_poles_4));
% z_poles_N = z_polesCorner_N*corners_N;

end