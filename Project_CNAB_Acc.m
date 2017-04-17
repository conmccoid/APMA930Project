%% Project - Burgers Equation with CNAB

%%
%===Chebyshev Collocation, CNAB===%

%---Parameters---%
v = 0.5;
T = 1;

%---Discretization (changing)---%
dti = logspace(-2,-4,5);
Nj = 2:8; Nj = 2.^Nj;

ErrorCC = zeros(length(dti),length(Nj));
for i = 1:length(dti)
    dt = dti(i);
    for j = 1:length(Nj)
        N = Nj(j);

%---Matrices for CNAB---%
[Dm,xc] = chebdifmat(N,2,1);
I = eye(N+1);

CN2 = I - v*dt/2*Dm(:,:,2);
CN2(1,:) = I(1,:); CN2(end,:) = I(end,:);

CN1 = I + v*dt/2*Dm(:,:,2);
CN1(1,:) = zeros(1,N+1); CN1(end,:) = zeros(1,N+1);

AB = -dt/2*Dm(:,:,1);
AB(1,:) = zeros(1,N+1); AB(end,:) = zeros(1,N+1);

%---Initialization---%
IC = exactSoln(xc,0,v);
u_old = IC;
u = u_old + dt*v*Dm(:,:,2)*u_old - dt*u_old.*(Dm(:,:,1)*u_old);
t = dt;

%---Algorithm---%
while t<T
    t = t+dt;
    rhs = CN1*u + 3*u.*(AB*u) - u_old.*(AB*u_old);
    rhs(1) = exactSoln(xc(1),t,v);
    rhs(end) = exactSoln(xc(end),t,v);
    u_new = CN2 \ rhs;
    u_old = u;
    u = u_new;
    
    Exact = exactSoln(xc,t,v);
end

ErrorCC(i,j) = norm(Exact - u);

    end
end

% figure(1)
% subplot(1,2,1)
% loglog(Nj,ErrorCC(1,:),'*--',Nj,ErrorCC(2,:),'*--',Nj,ErrorCC(3,:),'*--',...
%     Nj,ErrorCC(4,:),'*--',Nj,ErrorCC(5,:),'*--','MarkerSize',10,'Linewidth',2)
% xlabel('Number of points')
% ylabel('2-norm of error')
% title('Chebyshev collocation, v = 0.01')
% legend(['\Delta t = ', num2str(dti(1))],['\Delta t = ', num2str(dti(2))],...
%     ['\Delta t = ', num2str(dti(3))],['\Delta t = ', num2str(dti(4))],...
%     ['\Delta t = ', num2str(dti(5))])
% set(gca,'FontSize',16,'Linewidth',2)
% 
% figure(2)
% subplot(1,2,1)
% loglog(dti,ErrorCC(:,end-2),'*--','MarkerSize',10,'Linewidth',2)
% xlabel('\Delta t')
% ylabel('2-norm of error')
% title('Chebyshev collocation, v = 0.01')
% set(gca,'FontSize',16,'Linewidth',2)
%%
%===Chebyshev Tau, CNAB===%

%---Parameters---%
v = 0.5;
T = 1;

%---Discretization (changing)---%
dti = logspace(-2,-4,5);
Nj = 2:8; Nj = 2.^Nj;

ErrorCT = zeros(length(dti),length(Nj));
for m = 1:length(dti)
    dt = dti(m);
    for n = 1:length(Nj)
        N = Nj(n);

%---Grid---%
[~,xc] = chebdifmat(N,2,1);

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

%---Matrices for CNAB---%
I = eye(N+1);
CN2 = I - v*dt/2*Dspec*Dspec;
CN2(end-1,:) = ones(1,N+1);
CN2(end,:) = (-ones(1,N+1)).^(0:N);

CN1 = I + v*dt/2*Dspec*Dspec;

%---Initialization---%
IC = exactSoln(xc,0,v);
u_old = DCT*IC;
f_old = -DCT*((iDCT*u_old).*(iDCT*Dspec*u_old));
u = u_old + dt*v*Dspec*Dspec*u_old + dt*f_old;
t = dt;

%---Algorithm---%
while t<T
    t = t+dt;
    f = -DCT*((iDCT*u).*(iDCT*Dspec*u));
    rhs = CN1*u + dt/2*(3*f - f_old);
    rhs(end-1) = exactSoln(xc(1),t,v);
    rhs(end) = exactSoln(xc(end),t,v);
    u_new = CN2 \ rhs;
    u = u_new;
    f_old = f;
    
    Exact = exactSoln(xc,t,v);
end

ErrorCT(m,n) = norm(Exact - iDCT*u);

    end
end

% figure(1)
% subplot(1,2,2)
% loglog(Nj,ErrorCT(1,:),'*--',Nj,ErrorCT(2,:),'*--',Nj,ErrorCT(3,:),'*--',...
%     Nj,ErrorCT(4,:),'*--',Nj,ErrorCT(5,:),'*--','MarkerSize',10,'Linewidth',2)
% xlabel('Number of points')
% ylabel('2-norm of error')
% title('Chebyshev tau, v = 0.01')
% legend(['\Delta t = ', num2str(dti(1))],['\Delta t = ', num2str(dti(2))],...
%     ['\Delta t = ', num2str(dti(3))],['\Delta t = ', num2str(dti(4))],...
%     ['\Delta t = ', num2str(dti(5))])
% set(gca,'FontSize',16,'Linewidth',2)
% 
% figure(2)
% subplot(1,2,2)
% loglog(dti,ErrorCT(:,end-2),'*--','MarkerSize',10,'Linewidth',2)
% xlabel('\Delta t')
% ylabel('2-norm of error')
% title('Chebyshev tau, v = 0.01')
% set(gca,'FontSize',16,'Linewidth',2)

%%
%===Fourier Galerkin, CNAB===%

%---Parameters---%
v = 0.5;
T = 1;

%---Discretization---%
dti = logspace(-2,-4,5);
Nj = 2:8; Nj = 2.^Nj;

ErrorFG = zeros(length(dti),length(Nj));
for m = 1:length(dti)
    dt = dti(m);
    for n = 1:length(Nj)
        N = Nj(n);

%---Grid---%
dx = 2*pi/N;
x = 0:dx:2*pi-dx; x = x';

%---Matrices---%
if mod(N,2)==0
    k = [0:ceil(N/2)-1 0 -floor(N/2)+1:-1];
else
    k = [0:ceil(N/2)-1 -floor(N/2):-1];
end
% k = 2*pi*k;

I = eye(N);
CN2 = I + dt*v/2*diag(k.^2);
CN1 = I - dt*v/2*diag(k.^2);

AB = 1i*diag(k);

%---Initialization---%
IC = exactPeriodicSoln(x+4*pi,0,v,100);
u_old = fft(IC);
f_old = fft( (ifft( u_old ) ).*(ifft( AB*u_old ) ) );
u = u_old - dt*v*diag(k.^2)*u_old - dt*f_old;
t = dt;

%---Algorithm---%
while t<T
    t = t+dt;
    f = fft( (ifft( u ) ).*(ifft( AB*u ) ) );
    u_new = CN2 \ (CN1*u - dt/2*(3*f - f_old));
    u_old = u;
    u = u_new;
    f_old = f;
    
    u_phys = ifft(u);
    u_phys = [u_phys ; u_phys(1)];
    
    Exact = exactPeriodicSoln([x;2*pi]+4*pi,t,v,100);
end

ErrorFG(m,n) = norm(u_phys - Exact);

    end
end

% figure(3)
% loglog(Nj,ErrorFG(1,:),'*--',Nj,ErrorFG(2,:),'*--',Nj,ErrorFG(3,:),'*--',...
%     Nj,ErrorFG(4,:),'*--',Nj,ErrorFG(5,:),'*--','MarkerSize',10,'Linewidth',2)
% xlabel('Number of points')
% ylabel('2-norm of error')
% title('Fourier Galerkin, v = 0.5')
% legend(['\Delta t = ', num2str(dti(1))],['\Delta t = ', num2str(dti(2))],...
%     ['\Delta t = ', num2str(dti(3))],['\Delta t = ', num2str(dti(4))],...
%     ['\Delta t = ', num2str(dti(5))])
% set(gca,'FontSize',16,'Linewidth',2)
% 
% figure(4)
% loglog(dti,ErrorFG(:,end-2),'*--','MarkerSize',10,'Linewidth',2)
% xlabel('\Delta t')
% ylabel('2-norm of error')
% title('Fourier Galerkin, v = 0.5')
% set(gca,'FontSize',16,'Linewidth',2)

%%
%===Error graphs===%

%---Space discretization---%
ind = length(dti);
figure(1)
loglog(Nj,ErrorCC(ind,:),'*--',Nj,ErrorCT(ind,:),'*--',Nj,ErrorFG(ind,:),'*--','Markersize',10,'Linewidth',2)
xlabel('Number of points')
ylabel('2-norm of error')
title('\Deltat = 1e-4')
legend('Chebyshev collocation','Chebyshev tau','Fourier Galerkin')
set(gca,'Fontsize',16,'linewidth',2)

%---Time discretization---%
ind = length(Nj);
figure(2)
loglog(dti,ErrorCC(:,ind),'*--',dti,ErrorCT(:,ind),'*--',dti,ErrorFG(:,ind),'*--','Markersize',10,'Linewidth',2)
xlabel('Number of points')
ylabel('2-norm of error')
title('N = 256')
legend('Chebyshev collocation','Chebyshev tau','Fourier Galerkin')
set(gca,'Fontsize',16,'linewidth',2)

%---Linear fits---%
fitobject = fit(log(dti'),log(ErrorCC(:,ind)),'poly1');
orderCC = fitobject.p1;
fitobject = fit(log(dti'),log(ErrorCT(:,ind)),'poly1');
orderCT = fitobject.p1;
fitobject = fit(log(dti'),log(ErrorFG(:,ind)),'poly1');
orderFG = fitobject.p1;

disp(['Chebyshev collocation: ',num2str(orderCC),' Chebyshev tau: ',num2str(orderCT),' Fourier Galerkin: ',num2str(orderFG)])