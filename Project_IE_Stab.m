%% Project - Stability for IMEX Euler
% Finds stability region for linearized IMEX Euler methods

%%
%===CC Stability===%
alpha = logspace(0,-8,100);
% beta = logspace(0,-4,100);
beta = -1:1e-2:1;
stabCC = zeros(length(beta),length(alpha));

%---Grid---%
N = 64;
[Dm,~] = chebdifmat(N,2,1);
I = eye(N+1);

for i = 1:length(alpha)
    for j = 1:length(beta)
        A = I - alpha(i)*Dm(:,:,2);
        A(1,:) = I(1,:); A(end,:) = I(end,:);

        B = I - beta(j)*Dm(:,:,1);
        B(1,:) = zeros(1,N+1); B(end,:) = zeros(1,N+1);
        
        stabCC(j,i) = max(abs(eig(inv(A)*B)));
    end
end

%---Visualize---%
figure
pcolor(alpha,beta,stabCC)
shading flat
colormap jet
colorbar
caxis([min(min(stabCC)) max(max(stabCC))])
hold on
contour(alpha,beta,stabCC,[1 1],'k','Linewidth',10)
xlabel('v*dt')
ylabel('dt*c')
title('Chebyshev Collocation')
set(gca,'xscale','log','yscale','log','FontSize',26,'Linewidth',2)
hold off

%%
%===CT Stability===%
alpha = logspace(0,-8,100);
% beta = logspace(0,-4,100);
beta = -1:1e-2:1;
stabCT = zeros(length(beta),length(alpha));

%---Grid---%
N = 64;

%---Discrete Chebyshev Transform---%
i = 0:N; c = ones(size(i)); c(1) = 2; c(end) = 2;
DCT = bsxfun(@times,c',c);
DCT = (2/N)./DCT;
iDCT = cos(i'*i*pi/N);
DCT = DCT.*iDCT;

%---Spectral Differentiation---%
Dspec = zeros(N+1,N+1);
for i = 1:N+1
    for j = i+1:2:N+1
        if i==1
            Dspec(i,j) = (j-1);
        else
            Dspec(i,j) = 2*(j-1);
        end
    end
end

%---Matrices---%
I = eye(N+1);

for i = 1:length(alpha)
    for j = 1:length(beta)
        A = I - alpha(i)*Dspec*Dspec;
        A(end-1,:) = ones(1,N+1); A(end,:) = (-ones(1,N+1)).^(0:N);

        B = I - beta(j)*Dspec;
        B(end-1,:) = zeros(1,N+1); B(end,:) = zeros(1,N+1);
        
        stabCT(j,i) = max(abs(eig(inv(A)*B)));
    end
end

%---Visualize---%
figure
pcolor(alpha,beta,stabCT)
shading flat
colormap jet
colorbar
caxis([min(min(stabCT)) max(max(stabCT))])
hold on
contour(alpha,beta,stabCT,[1 1],'k','Linewidth',10)
xlabel('v*dt')
ylabel('dt*c')
title('Chebyshev Tau')
set(gca,'xscale','log','yscale','log','FontSize',26,'Linewidth',2)
hold off

%%
%===FG Stability===%
alpha = logspace(0,-8,100);
% beta = logspace(0,-4,100);
beta = -1:1e-2:1;
stabFG = zeros(length(beta),length(alpha));

%---Grid---%
N = 65;

%---Matrices---%
if mod(N,2)==0
    k = [0:ceil(N/2)-1 0 -floor(N/2)+1:-1];
else
    k = [0:ceil(N/2)-1 -floor(N/2):-1];
end
k = 2*pi*k;

I = eye(N);

for i = 1:length(alpha)
    for j = 1:length(beta)
        A = I + alpha(i)*diag(k.^2);
        B = I - beta(j)*1i*diag(k);
        
        lambda = eig(inv(A)*B);
%         plot(real(lambda),imag(lambda),'r.','MarkerSize',10)
%         hold on
%         theta = -1:0.0001:1;
%         plot(cos(pi*theta),sin(pi*theta),'k--','linewidth',2)
%         hold off
%         axis equal
%         xlabel('Real part of eigenvalues')
%         ylabel('Imaginary part')
%         pause(0.01)
        
        stabFG(j,i) = max(abs(lambda));
    end
end

%---Visualize---%
figure
pcolor(alpha,beta,stabFG)
shading flat
colormap jet
colorbar
caxis([min(min(stabFG)) max(max(stabFG))])
hold on
contour(alpha,beta,stabFG,[1+eps 1+eps],'k','Linewidth',10)
xlabel('v*dt')
ylabel('c*dt')
title('Fourier Galerkin')
set(gca,'xscale','log','yscale','log','FontSize',26,'Linewidth',2)
hold off

%%
%===Combined Visualization===%
figure

%---Chebyshev collocation---%
subtightplot(3,1,1)
pcolor(alpha,beta,stabCC)
shading flat
colormap jet
colorbar
caxis([min(min(stabCC)) max(max(stabCC))])
hold on
contour(alpha,beta,stabCC,[1 1],'k','Linewidth',10)
ylabel('c\Deltat')
title('v\Deltat','Fontweight','normal')
set(gca,'xscale','log','XTick',[],'FontSize',20,'Linewidth',2)
hold off

%---Chebyshev tau---%
subtightplot(3,1,2)
pcolor(alpha,beta,stabCT)
shading flat
colormap jet
colorbar
caxis([min(min(stabCT)) max(max(stabCT))])
hold on
contour(alpha,beta,stabCT,[1 1],'k','Linewidth',10)
ylabel('Chebyshev Tau')
set(gca,'xscale','log','XTick',[],'YTick',[],'FontSize',20,'Linewidth',2)
hold off

%---Fourier Galerkin---%
subtightplot(3,1,3)
pcolor(alpha,beta,stabFG)
shading flat
colormap jet
colorbar
caxis([min(min(stabFG)) max(max(stabFG))])
hold on
contour(alpha,beta,stabFG,[1+eps 1+eps],'k','Linewidth',10)
ylabel('Fourier Galerkin')
set(gca,'xscale','log','XTick',[1e-7 1e-5 1e-3 1e-1],'YTick',[],'FontSize',20,'Linewidth',2)
hold off