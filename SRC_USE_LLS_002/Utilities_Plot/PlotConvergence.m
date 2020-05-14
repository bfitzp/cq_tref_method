function PlotConvergence(figPos,N_vec,bdyErr,A_cond,x_norm)
sqrt_N = sqrt(N_vec);

hFig = figure;
set(hFig, 'Position', figPos)
subplot(3,1,1)
hold on, grid on
plot(sqrt_N,bdyErr,'b');
plot(sqrt_N,bdyErr,'b*');
set(gca,'yscale','log')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('$||Ax-b||_\infty$','Interpreter','latex')

subplot(3,1,2)
hold on, grid on
plot(sqrt_N,A_cond,'b');
plot(sqrt_N,A_cond,'b*');
set(gca,'yscale','log')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('$cond(A)$','Interpreter','latex')

subplot(3,1,3)
hold on, grid on
plot(sqrt_N,x_norm,'b');
plot(sqrt_N,x_norm,'b*');
set(gca,'yscale','log')
xlabel('$\sqrt{N}$','Interpreter','latex')
ylabel('$||x||_2$','Interpreter','latex')
end