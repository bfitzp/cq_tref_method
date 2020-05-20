function [Z,G,tt] = GenSamplePointsAndBdyData(...
    w,wc,wr,wi,ww,dw,nw,pt,outward,g,polmax,scl,k_Nyq)

Z = [];           % col vector of sample points on boundary
G = [];           % col vector of boundary values at these points
tt = cell(nw,1);  % cell array of distances of sample points along each side

inpolygonc = @(z,w) inpolygon(real(z), ...
    imag(z),real(w),imag(w));

%----------------------------------------------------------------------
% Generate sets of distances along each edge
for k = 1:nw
    nk = polmax;                                  % no. of poles at this corner
    sk = sqrt(1:nk) - sqrt(nk);
    dk = exp(2*sk); dk = scl*dk;                  % stronger clustering near corner
    dk = dk(dk>1e-15*scl);                        % remove poles too close to corner
    polk = w(k) + outward(k)*dk;                  % poles near this corner
    ii = find(inpolygonc(polk(dk>1e-12*scl),ww),1); % work around inaccuracy
    if length(ii)>0                               % don't allow poles in Omega
        dk = dk(1:ii-2); polk = polk(1:ii-2);
    end

    dvec = [(1/3)*dk (2/3)*dk dk];                % finer pts for bndry sampling
    
%     tt_N = 5;
%     tt{k} = [tt{k} dvec(dvec<dw(k)) ...  
%         linspace(0,dw(k),max(tt_N,nk))];   
%     j = mod(k-2,nw)+1;                            % index of last corner
%     tt{j} = [tt{j} dw(j)-dvec(dvec<dw(j))];       % likewise in other direction

    tt_ = 2*pi*(0:k_Nyq)/k_Nyq;
    tt_ = tt_(tt_<dw(k));
    tt{k} = [tt{k} dvec(dvec<dw(k)) tt_];   
    j = mod(k-2,nw)+1;
    tt{j} = [tt{j} dw(j)-dvec(dvec<dw(j))];       % likewise in other direction
end

%----------------------------------------------------------------------
% Generate sample points and boundary data
for k = 1:nw
    tt{k} = sort(tt{k}(:));
    tk = tt{k}; pk = pt{k};                       % abbrevations
    Z = [Z; pk(tk)];                              % sample pts on side k
    G = [G; g{k}(pk(tk))];                        % boundary data at these pts
%     Z = [Z; pk];                              % sample pts on side k
%     G = [G; g{k}(pk)];                        % boundary data at these pts
    h = 1e-4;                                     % 4-pt trapezoidal rule
end
 
% LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
% fs = 9; PO = 'position'; FW = 'fontweight'; NO = 'normal';
% sx = linspace(wr(1),wr(2),200); sy = linspace(wi(1),wi(2),200);
% [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;
% ax = [wr(1:2); wi(1:2)] + .2*scl*[-1 1 -1 1]';
% axwide = [wr(1:2); wi(1:2)] + 1.1*scl*[-1 1 -1 1]';
%     
%     
% figure
% hold on, grid on
% plot(ww,'k-',LW,1)
% scatter(real(Z),imag(Z),'r+',LW,1)
% set(gca,FS,fs-1)
% plot(real(wc),imag(wc),'k.',MS,6)
% axis(ax)
% axis equal

end