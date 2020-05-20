%----------------------------------------------------------------------
% Hardcoded for the moment...
%----------------------------------------------------------------------
function [z_poles,N_1]...
    = GenPoles_New(...
    z_polesCorner_N,z_polesOffset,sigma,...
    w, cornerInds)

corners_N = length(w);

j = (0:z_polesCorner_N-1).';
polesPos = exp(-sigma*j/sqrt(z_polesCorner_N));
polesPos = polesPos - polesPos(end);

poles_MinDist = abs(polesPos(end) - polesPos(end-1));

corners_N = length(cornerInds);

z_poles = zeros(z_polesCorner_N,corners_N);

for i=1:corners_N   
%     if i == 1
%         a_2 = w(cornerInds(2))-w(cornerInds(1));
%         a_1 = w(cornerInds(end))-w(cornerInds(1));
%     elseif i < corners_N
%         a_2 = w(cornerInds(i+1))-w(cornerInds(i));
%         a_1 = w(cornerInds(i-1))-w(cornerInds(i));
%     else
%         a_2 = w(cornerInds(1))-w(cornerInds(end));
%         a_1 = w(cornerInds(end-1))-w(cornerInds(end));
%     end
    if i == 1
        a_2 = w(2)-w(1);
        a_1 = w(end)-w(1);
    elseif i < corners_N
        a_2 = w(i+1)-w(i);
        a_1 = w(i-1)-w(i);
    else
        a_2 = w(1)-w(end);
        a_1 = w(end-1)-w(end);
    end
    
    c = norm(a_2)*a_1 + norm(a_1)*a_2;
    c_nrmlzd = c/norm(c);
    
    %----------------------------------------------------------------------
    % Place the poles in the appropriate direction
    cp = cross([real(a_2),imag(a_2),0],[real(a_1),imag(a_1),0]);    
    if cp >= 0
        poleSign = +1;
    else
        poleSign = -1;
    end
    z_poles(:,i) = w(i) + poleSign*(...
        polesPos*c_nrmlzd + poles_MinDist*c_nrmlzd); % + z_polesOffset*c_nrmlzd);
end

z_poles = z_poles(:);

N_1 = corners_N*z_polesCorner_N;

end   
    

% fprintf('i: %d /// (a_2,a_1): (%.2f,%.2f) /// cp: %.2f \n',i,a_2,a_1,cp(3))
% quiver(real(w(i)), imag(w(i)), real(a_2), imag(a_2), 'b')
% quiver(real(w(i)), imag(w(i)), real(a_1), imag(a_1), 'r')
% scatter(real(z_poles(:,i)),imag(z_poles(:,i)))












% 
% %----------------------------------------------------------------------
% % Hardcoded for the moment...
% %----------------------------------------------------------------------
% function [z_poles,N_1]...
%     = GenPoles_New(...
%     z_polesCorner_N,z_polesOffset,sigma,...
%     w, cornerInds)
% 
% corners_N = length(w);
% 
% j = (0:z_polesCorner_N-1).';
% 
% w_N = length(w);
% z_poles = zeros(z_polesCorner_N,w_N);
% for i=1:w_N
%     if i == 1
%         a_2 = w(2)-w(1);
%         a_1 = w(end)-w(1);
%     elseif i < w_N
%         a_2 = w(i+1)-w(i);
%         a_1 = w(i-1)-w(i);
%     else
%         a_2 = w(1)-w(end);
%         a_1 = w(end-1)-w(end);
%     end
%     
%     c = norm(a_2)*a_1 + norm(a_1)*a_2;
%     c_nrmlzd = c/norm(c);
%     
%     %----------------------------------------------------------------------
%     % Place the poles in the appropriate direction
%     cp = cross([real(a_2),imag(a_2),0],[real(a_1),imag(a_1),0]);    
%     if cp >= 0
%         poleSign = +1;
%     else
%         poleSign = -1;
%     end
%     ppp = w(i) + poleSign*(exp(-sigma*j/sqrt(z_polesCorner_N))*c_nrmlzd + z_polesOffset*c_nrmlzd);
%     z_poles(:,i) = ppp;
% end
% 
% z_poles = z_poles(:);
% 
% N_1 = corners_N*z_polesCorner_N;
% 
% end   
%     
% 
% % fprintf('i: %d /// (a_2,a_1): (%.2f,%.2f) /// cp: %.2f \n',i,a_2,a_1,cp(3))
% % quiver(real(w(i)), imag(w(i)), real(a_2), imag(a_2), 'b')
% % quiver(real(w(i)), imag(w(i)), real(a_1), imag(a_1), 'r')
% % scatter(real(z_poles(:,i)),imag(z_poles(:,i)))