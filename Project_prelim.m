%%Project - Preliminary
%%
%===Heat Equation, Chebyshev Collocation===%

%---Grid---%
N = 64;
[Dm,xc] = chebdifmat(N,2,1);
T = 1;
dt = 9/N^2;

%---Matrices for Algorithms---%
I = eye(N+1);
FE = I + dt*Dm(:,:,2);
FE(1,:) = I(1,:); FE(end,:) = I(end,:);

BE = I - dt*Dm(:,:,2);
BE(1,:) = I(1,:); BE(end,:) = I(end,:);

AB2 = I + 3*dt*Dm(:,:,2)/2;

CN1 = I + dt*Dm(:,:,2)/2;
CN2 = I - dt*Dm(:,:,2)/2;
CN1(1,:) = I(1,:); CN1(end,:) = I(end,:);
CN2(1,:) = I(1,:); CN2(end,:) = I(end,:);

%---Initial conditions---%
t = 0;
IC = sin(pi*xc);
u_old = IC;
u = u_old + dt*Dm(:,:,2)*u_old;

%---Algorithm---%
while t<T
    t = t+dt;
    
    %---Forward Euler---%
%     u = FE*u;

    %---Backward Euler---%
%     u = BE\u;
    
    %---2nd Order Adams-Bashforth---%
%     u = AB2*u - Dm(:,:,2)*u_old/2;
%     u(1) = 0; u(end) = 0;
%     u_old = u;

    %---Crank-Nicolson---%
    u = CN2\(CN1*u);

    Exact = IC*exp(-t*pi^2);
    plot(xc,u,xc,Exact)
    axis([-1,1,-1,1])
    pause(0.01)
end

%%
%===-u_xx + (y^2)u = f, Chebyshev Tau===%

%---Grid---%
N = 64;
[~,xc] = chebdifmat(N,1,1);

I = eye(N+1);

%---Discrete Chebyshev Transform---%
i = 0:N; c = ones(size(i)); c(1) = 2; c(end) = 2;
DCT = bsxfun(@times,c',c);
DCT = (2/N)./DCT;
iDCT = cos(i'*i*pi/N);
DCT = DCT.*iDCT;

%---Spectral Differentiation---%
Dspec = zeros(N+1,N+1);
for i = 1:N
    for j = i+1:2:N+1
        if i==1
            Dspec(i,j) = (j-1);
        else
            Dspec(i,j) = 2*(j-1);
        end
    end
end

%---Problem Parameters---%
y = 1;
exact = sin(pi*xc);
f = (pi^2 + y^2)*sin(pi*xc);
F = DCT*f; F(end-1:end) = [0;0];
A = (y^2)*I - Dspec*Dspec;
A(end-1,:) = ones(1,N+1);
A(end,:) = (-ones(1,N+1)).^(0:N);

%---Initial Guess---%
% u = zeros(N+1,1);
% err = 1;
% tol = 1e-9;

%---Algorithm---%
u = A\F;

plot(xc,iDCT*u,'bo',xc,exact,'r*')