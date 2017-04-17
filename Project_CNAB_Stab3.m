%% Project - Stability for CNAB

%%
%===Chebyshev collocation===%

%---Grid---%
alpha = logspace(0,-4,100);
beta = logspace(0,-6,100);
% beta = -1:1e-2:1;
stabCC = zeros(length(beta),length(alpha));

%---Standard matrices---%
[Dm,~] = chebdifmat(64,2,1);
I = eye(65);

%---Algorithm matrices---%
for i = 1:length(alpha)
    for j = 1:length(beta)
        
        A = I - alpha(i)*beta(j)/2*Dm(:,:,2);
        A(1,:) = I(1,:); A(end,:) = I(end,:);
        A = inv(A);

        B = I + alpha(i)*beta(j)/2*Dm(:,:,2) - 3*beta(j)/2*Dm(:,:,1);
        B(1,:) = zeros(65,1); B(end,:) = zeros(65,1);

        C = beta(j)/2*Dm(:,:,1);
        C(1,:) = zeros(65,1); C(end,:) = zeros(65,1);

        Alpha = [ zeros(65) I ; A*C A*B ];

        stabCC(j,i) = max(abs(eig(Alpha)));

    end
end

%---Visualization---%
figure
pcolor(alpha,beta,stabCC)
shading flat
colormap jet
colorbar
caxis([min(min(stabCC)) max(max(stabCC))])
hold on
contour(alpha,beta,stabCC,[1 1],'k','Linewidth',10)
xlabel('v')
ylabel('\Deltat')
title('Chebyshev Collocation')
set(gca,'xscale','log','yscale','log','FontSize',26,'Linewidth',2)
hold off

%%
%===Chebyshev tau===%
alpha = logspace(0,-4,100);
beta = logspace(0,-6,100);
% beta = -1:1e-2:1;
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
        A = I - alpha(i)*beta(j)/2*Dspec*Dspec;
        A(end-1,:) = ones(1,N+1); A(end,:) = (-ones(1,N+1)).^(0:N);
        A = inv(A);
        
        B = I + alpha(i)*beta(j)/2*Dspec*Dspec - 3*beta(j)/2*Dspec;
        B(end-1,:) = zeros(1,N+1); B(end,:) = zeros(1,N+1);
        
        C = beta(j)/2*Dspec;
        C(end-1,:) = zeros(1,N+1); C(end,:) = zeros(1,N+1);
        
        Alpha = [ zeros(N+1) I ; A*C A*B ];
        
        stabCT(j,i) = max(abs(eig(Alpha)));
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
xlabel('v')
ylabel('\Deltat')
title('Chebyshev Tau')
set(gca,'xscale','log','yscale','log','FontSize',26,'Linewidth',2)
hold off

%%
%===FG Stability===%
alpha = logspace(0,-4,100);
beta = logspace(0,-6,100);
% beta = -1:1e-2:1;
stabFG = zeros(length(beta),length(alpha));

%---Grid---%
N = 65;

%---Matrices---%
if mod(N,2)==0
    k = [0:ceil(N/2)-1 0 -floor(N/2)+1:-1];
else
    k = [0:ceil(N/2)-1 -floor(N/2):-1];
end
% k = 2*pi*k;

I = eye(N);

for i = 1:length(alpha)
    for j = 1:length(beta)
        A = I + alpha(i)/2*diag(k.^2); A = inv(A);
        B = I - alpha(i)/2*diag(k.^2) - 3*beta(j)*1i/2*diag(k);
        C = beta(j)*1i/2*diag(k);
        
        Alpha = [ zeros(N) I ; A*C A*B ];
                
        stabFG(j,i) = max(abs(eig(Alpha)));
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
xlabel('v')
ylabel('\Deltat')
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
ylabel('\Deltat')
title('v','Fontweight','normal')
set(gca,'xscale','log','yscale','log','XTick',[],'FontSize',20,'Linewidth',2)
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
set(gca,'xscale','log','yscale','log','XTick',[],'YTick',[],'FontSize',20,'Linewidth',2)
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
set(gca,'xscale','log','yscale','log','XTick',[1e-3 1e-1],'YTick',[],'FontSize',20,'Linewidth',2)
hold off