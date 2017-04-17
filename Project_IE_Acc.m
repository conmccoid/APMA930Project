%% Project - Accuracy and Efficiency(?) for IMEX Euler
% Calculates accuracy and efficiency(?) for IMEX Euler methods

% Accuracy will likely change with N (space discretization), t (time), dt
% (time discretization), and v (physical parameter)
% Final time T and parameter v should remain fixed, so that we only check
% how discretization affects accuracy

%%
%===CC Accuracy===%

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

A = I - dt*v*Dm(:,:,2);
A(1,:) = I(1,:); A(end,:) = I(end,:);

% ---Initialization---%
IC = exactSoln(xc,0,v);
u = IC;
t = 0;

%---Algorithm---%
while t<T
    t = t+dt;
    rhs = u - dt*(u.*(Dm(:,:,1)*u));
    rhs(1) = exactSoln(xc(1),t,v);
    rhs(end) = exactSoln(xc(end),t,v);
    u = A \ rhs;
    
    Exact = exactSoln(xc,t,v);
end

ErrorCC(i,j) = norm(Exact - u);

    end
end

figure(1)
subplot(1,2,1)
loglog(Nj,ErrorCC(1,:),'*--',Nj,ErrorCC(2,:),'*--',Nj,ErrorCC(3,:),'*--',...
    Nj,ErrorCC(4,:),'*--',Nj,ErrorCC(5,:),'*--','MarkerSize',10,'Linewidth',2)
xlabel('Number of points')
ylabel('2-norm of error')
title('Chebyshev collocation, v = 0.01')
legend(['\Delta t = ', num2str(dti(1))],['\Delta t = ', num2str(dti(2))],...
    ['\Delta t = ', num2str(dti(3))],['\Delta t = ', num2str(dti(4))],...
    ['\Delta t = ', num2str(dti(5))])
set(gca,'FontSize',16,'Linewidth',2)

figure(2)
subplot(1,2,1)
loglog(dti,ErrorCC(:,end-2),'*--','MarkerSize',10,'Linewidth',2)
xlabel('\Delta t')
ylabel('2-norm of error')
title('Chebyshev collocation, v = 0.01')
set(gca,'FontSize',16,'Linewidth',2)

%%
%===CT Accuracy===%

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

%---Matrices---%
I = eye(N+1);

A = I - v*dt*Dspec*Dspec;
A(end-1,:) = ones(1,N+1);
A(end,:) = (-ones(1,N+1)).^(0:N);

%---Initialization---%
IC = exactSoln(xc,0,v);
u = DCT*IC;
t = 0;

%---Algorithm---%
while t<T
    t = t+dt;
    rhs = u - dt * ( DCT * ( (iDCT*u).*(iDCT*Dspec*u) ) );
    rhs(end-1) = exactSoln(xc(1),t,v);
    rhs(end) = exactSoln(xc(end),t,v);
    u = A \ rhs;
    
    Exact = exactSoln(xc,t,v);
end

ErrorCT(m,n) = norm(Exact - iDCT*u);

    end
end

figure(1)
subplot(1,2,2)
loglog(Nj,ErrorCT(1,:),'*--',Nj,ErrorCT(2,:),'*--',Nj,ErrorCT(3,:),'*--',...
    Nj,ErrorCT(4,:),'*--',Nj,ErrorCT(5,:),'*--','MarkerSize',10,'Linewidth',2)
xlabel('Number of points')
ylabel('2-norm of error')
title('Chebyshev tau, v = 0.01')
legend(['\Delta t = ', num2str(dti(1))],['\Delta t = ', num2str(dti(2))],...
    ['\Delta t = ', num2str(dti(3))],['\Delta t = ', num2str(dti(4))],...
    ['\Delta t = ', num2str(dti(5))])
set(gca,'FontSize',16,'Linewidth',2)

figure(2)
subplot(1,2,2)
loglog(dti,ErrorCT(:,end-2),'*--','MarkerSize',10,'Linewidth',2)
xlabel('\Delta t')
ylabel('2-norm of error')
title('Chebyshev tau, v = 0.01')
set(gca,'FontSize',16,'Linewidth',2)

%%
%===Fourier Galerkin===%

%---Parameters---%
v = 0.01;
T = 1;

%---Discretization---%
dti = logspace(-2,-4,5);
Nj = 2:8; Nj = 2.^Nj;

ErrorFG = zeros(length(dti),length(Nj));
for m = 1:length(dti)
    dt = dti(m);
    for n = 1:length(Nj)
        N = Nj(n);

% N = 64;
% dt = 0.01;

%---Grid---%
dx = 2*pi/N;
x = 0:dx:2*pi-dx; x = x';

%---Matrices---%
if mod(N,2)==0
    k = [0:ceil(N/2)-1 0 -floor(N/2)+1:-1];
else
    k = [0:ceil(N/2)-1 -floor(N/2):-1];
end

I = eye(N);
A = I + dt*v*diag(k.^2);
B = 1i*diag(k);

%---Initialization---%
IC = exactPeriodicSoln(x+4*pi,0,v,100);
u = fft(IC);
t = 0;

%---Algorithm---%
while t<T
    t = t+dt;
    rhs = u - dt*fft( (ifft(u).*ifft(B*u)) );
    u = A \ rhs;
    
    u_phys = ifft(u,'symmetric');
    u_phys = [u_phys ; u_phys(1)];
    
    Exact = exactPeriodicSoln([x+4*pi;6*pi],t,v,100);
    
%     figure(3)
%     plot([x;2*pi], u_phys,'r--',[x;2*pi], exactPeriodicSoln([x+4*pi;6*pi],t,v,100),'b')
%     axis([0,2*pi,-3,3])
%     pause(0.01)
end

ErrorFG(m,n) = norm(u_phys - Exact);

    end
end

figure(3)
loglog(Nj,ErrorFG(1,:),'*--',Nj,ErrorFG(2,:),'*--',Nj,ErrorFG(3,:),'*--',...
    Nj,ErrorFG(4,:),'*--',Nj,ErrorFG(5,:),'*--','MarkerSize',10,'Linewidth',2)
xlabel('Number of points')
ylabel('2-norm of error')
title('Fourier Galerkin, v = 0.01')
legend(['\Delta t = ', num2str(dti(1))],['\Delta t = ', num2str(dti(2))],...
    ['\Delta t = ', num2str(dti(3))],['\Delta t = ', num2str(dti(4))],...
    ['\Delta t = ', num2str(dti(5))])
set(gca,'FontSize',16,'Linewidth',2)

figure(4)
loglog(dti,ErrorFG(:,end-2),'*--','MarkerSize',10,'Linewidth',2)
xlabel('\Delta t')
ylabel('2-norm of error')
title('Fourier Galerkin, v = 0.01')
set(gca,'FontSize',16,'Linewidth',2)