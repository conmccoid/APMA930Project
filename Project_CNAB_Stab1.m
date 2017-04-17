%% Project - Stability for CNAB

%%
%===Chebyshev collocation===%

%---Parameters---%
v = 0.01;
c = 1;
dt = 0.01;

%---Eigenvalues---%
[Dm,~] = chebdifmat(64,1,1);
% lambdaCC = eig(Dm(2:end-1,2:end-1,1));
lambdaCC = eig(Dm);

%---Grid---%
[lambda_real,lambda_imag] = meshgrid(-100:0.5:100,-100:0.5:100);
lambda = lambda_real + 1i*lambda_imag;
a = 1 - v*dt/2*lambda.^2;
b = 1 + v*dt/2*lambda.^2 - 3*c*dt/2*lambda;
d = c*dt/2*lambda;

% [m,n] = size(lambda);
% stabCC = zeros(m,n);
% for i = 1:m
%     for j = 1:n
%         stabCC(i,j) = max(abs(roots([a(i,j),-b(i,j),-d(i,j)])));
%     end
% end

r1 = ( b + sqrt(b.^2 + 4*a.*d) )./(2*a);
r2 = ( b - sqrt(b.^2 + 4*a.*d) )./(2*a);
stabCC = max(abs(r1),abs(r2));

%---Visualization---%
figure(1)
pcolor(lambda_real,lambda_imag,stabCC)
shading flat
colormap jet
colorbar
caxis([min(min(stabCC)) max(max(stabCC))])
hold on
contour(lambda_real,lambda_imag,stabCC,[1 1],'k','Linewidth',10)
plot(real(lambdaCC),imag(lambdaCC),'m.','MarkerSize',25)
% plot(real(lambdaCT),imag(lambdaCT),'k^','MarkerSize',10)
xlabel('Real')
ylabel('Imag')
set(gca,'FontSize',26,'Linewidth',10)
hold off