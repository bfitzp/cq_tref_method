function PlotCQFreqs(f,T,M,p)
t = linspace(0,T,M+1); % t = linspace(0,T,k*M+1);
F = f(t);
dt = T/M; % kappaCQ=T/(k*M);

N  = size(F,2)-1;  
R = eps^(0.5/(N+1));
omega = exp(2*pi*1i/(N+1));

figure
hold on, grid on
for l=0:floor((N+1)/2)
    k_l = 1i*p(R*omega^(-l))/dt;
    scatter(real(k_l), imag(k_l),'b')
end
axis equal

end


% doPlotFreqs = 0;
% if doPlotFreqs
%     figure
%     hold on
%     for l=0:N_F
%     % for l=0:floor((N_F+1)/2)
%         k_l = 1i*p(R*omega^(-l))/kappaCQ;
%
%         if l<=floor((N_F+1)/2)
%             dotColor = [0.66,0.20,0.18];
%         else
%             dotColor = [0,0,0];
%         end
%         scatter(real(k_l), imag(k_l),1,'MarkerEdgeColor',dotColor,...
%                   'MarkerFaceColor',dotColor,...
%                   'LineWidth',0.1)
%     end
%     xlabel('real$(k_l)$','Interpreter','latex')
%     ylabel('imag$(k_l)$','Interpreter','latex')
%     grid on
%     set(gca,'Color',[.98,.98,.98])
%     set(gcf,'Color',[.98,.98,.98])
%     set(gca,'FontSize',16)
%     PrintFigure('complex_wavenumbers',6,4)
% end