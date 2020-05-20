%----------------------------------------------------------------------
% Hardcoded for the moment...
%----------------------------------------------------------------------
function [z_samp,z_samp_weights,cornerInds,M] = GenSamplePoints_New(...
    w,edgePanels_N,refine_N)

inc = ((0:edgePanels_N-1)/edgePanels_N).';
inc_N = length(inc);

%----------------------------------------------------------------------
% Create sample points
w_N = length(w);
v = zeros(w_N*inc_N,1);
for j=1:w_N-1
    edgeDir = w(j+1)-w(j);
    v((j-1)*inc_N+1:j*inc_N) = w(j) + edgeDir*inc;
end

edgeDir = w(1)-w(w_N);
v((w_N-1)*inc_N+1:end) = w(end) + edgeDir*inc;

%----------------------------------------------------------------------
% Set the new corner indices
cornerInds = (1 + (0:w_N-1)*edgePanels_N).';

%----------------------------------------------------------------------
% Refine the panels closest to the corners
for j=1:refine_N
    [v, cornerInds] = refineGrid(v,w_N,edgePanels_N,cornerInds);
end

z_samp = v;
M = length(v);

z_samp_weights = abs(diff(z_samp));
z_samp_weights = [z_samp_weights; abs(z_samp(1)-z_samp(end))];
z_samp_weights = sqrt(z_samp_weights);
z_samp_weights = z_samp_weights.^0;
end

function [v_new, cornerInds_new] = refineGrid(v,w_N,edgePanels_N,cornerInds)
v_new = [];
for j=1:w_N+1
    if j==1
        refDir_1_b = v(2)-v(1);
        pNew_1_b = v(1) + refDir_1_b/2;
        
        v_new = [v(1); pNew_1_b];
    elseif j < w_N+1
        ind_prev = cornerInds(j-1);
        ind = cornerInds(j);
        refDir_a = v(ind-1)-v(ind);
        pNew_a = v(ind) + refDir_a/2;

        refDir_b = v(ind+1)-v(ind);
        pNew_b = v(ind) + refDir_b/2;
        
        v_new = [v_new; v(ind_prev+1:ind-1); pNew_a; v(ind); pNew_b];
    else
        ind_prev = cornerInds(j-1);
        ind = 1 + (j-1)*edgePanels_N;
        
        refDir_1_a = v(end)-v(1);
        pNew_1_a = v(1) + refDir_1_a/2;

        v_new = [v_new; v(ind_prev+1:end); pNew_1_a];
    end
end

cornerInds_new = cornerInds + (2*(0:w_N-1)).';

end










% figure
% hold on, grid on
% scatter(real(z_samp), imag(z_samp), 'b')





% 
% 
% clf
% figure
% hold on, grid on
% scatter(real(v_1), imag(v_1), 'b+')
% 
% [v_2, cornerInds_2] = refineGrid(v_1,w_N,edgePanels_N,cornerInds_1);
% 
% clf
% figure
% hold on, grid on
% scatter(real(v_2), imag(v_2), 'b+')
% 
% [v_3, cornerInds_3] = refineGrid(v_2,w_N,edgePanels_N,cornerInds_2);
% 
% clf
% figure
% hold on, grid on
% scatter(real(v_3), imag(v_3), 'b+')










% figure
% hold on, grid on
% scatter(real(v), imag(v), 'b')
% scatter(real(v_new), imag(v_new), 'r+')
% abc = 1;


% refDir_2_a = v(3)-v(4);
% pNew_2_a = v(4) + refDir_2_a/2;
% 
% refDir_2_b = v(5)-v(4);
% pNew_2_b = v(4) + refDir_2_b/2;
% 
% 
% refDir_3_a = v(6)-v(7);
% pNew_3_a = v(7) + refDir_3_a/2;
% 
% refDir_3_b = v(8)-v(7);
% pNew_3_b = v(7) + refDir_3_b/2;
% 
% 
% refDir_4_a = v(9)-v(10);
% pNew_4_a = v(10) + refDir_4_a/2;
% 
% refDir_4_b = v(11)-v(10);
% pNew_4_b = v(10) + refDir_4_b/2;


% 
% 
% v = [...
% % j=1                           v(1);  pNew_1_b;
% % j=2       v(2:3);   pNew_2_a; v(4);  pNew_2_b;
% % j=3       v(5:6);   pNew_3_a; v(7);  pNew_3_b;...
% % j=4       v(8:9);   pNew_4_a; v(10); pNew_4_b;...
% % j=5       v(11:12); pNew_1_a;
%     ];
% 
% cornerInds = (1 + (0:w_N-1)*(edgePanels_N+2)).';
% 
% % 1       1
% % 4       6
% % 7       11
% % 10      16
% 
% % cornerInds(2) = cornerInds(2) + 1;
% 
% 
% % refDir = v(5)-v(4);
% % pNew = v(4) + refDir/2;
% % v = [v(1:4); pNew; v(5:end)];
% % 
% 
% figure
% scatter(real(v), imag(v), 'b')
% 
% 
% aaa = 1;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% incFrac = (1:z_samp_N)/z_samp_N;
% 
% %----------------------------------------------------------------------
% % Sample points - w_1
% z_samp_1 = GenerateSamplePoints_(...
%     z_poles(:,1),w(1),incFrac,...
%     z_polesCorner_N,z_samp_N,1,1);
% 
% %----------------------------------------------------------------------
% % Sample points - w_2
% z_samp_2 = GenerateSamplePoints_(...
%     z_poles(:,2),w(2),incFrac,...
%     z_polesCorner_N,z_samp_N,1,-1);
% 
% %----------------------------------------------------------------------
% % Sample points - w_3
% z_samp_3 = GenerateSamplePoints_(...
%     z_poles(:,3),w(3),incFrac,...
%     z_polesCorner_N,z_samp_N,-1,-1);
% 
% %----------------------------------------------------------------------
% % Sample points - w_4
% z_samp_4 = GenerateSamplePoints_(...
%     z_poles(:,4),w(4),incFrac,...
%     z_polesCorner_N,z_samp_N,-1,1);
% 
% z_samp = [...
%     z_samp_1;...
%     z_samp_2;...
%     z_samp_3;...
%     z_samp_4;...
%     ];
% 
% z_samp_N = length(z_samp);
% M = z_samp_N;
% 
% end
% 
% function z_samp = GenerateSamplePoints_(...
%     z_poles_vec,w,incFrac,...
%     z_polesCorner_N,z_samp_N,a,b)
% delta = abs(z_poles_vec - w);
% inc = delta*incFrac;
% z_samp_h = w*ones(z_polesCorner_N,z_samp_N) + 4*a*inc;
% z_samp_v = w*ones(z_polesCorner_N,z_samp_N) + 4*b*1i*inc;
% 
% z_samp = [z_samp_h(:);z_samp_v(:)];
% 
% threshold = 0.99;
% z_samp = z_samp(find(abs(z_samp-w) < threshold));
% 
% end