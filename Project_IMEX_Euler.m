%% Project - IMEX Euler

%%
%===Chebyshev Collocation===%

%---Parameters---%
v = 4e-4;

%---Grid---%
N = 64;
[Dm,xc] = chebdifmat(N,2,1);
T = 4;
dt = 1e-6;

%---Matrices for IMEX Euler---%
I = eye(N+1);

A = I - dt*v*Dm(:,:,2);
A(1,:) = I(1,:); A(end,:) = I(end,:);

% ---Initialization---%
% IC = sin(2*pi*xc) + sin(pi*xc)/2;
IC = sin(pi*xc);
% IC = exactSoln(xc,0,v);
u = IC;
t = dt;

%---Algorithm---%
while t<T
    t = t+dt;
    rhs = u - dt*(u.*(Dm(:,:,1)*u));
    rhs(1) = 0;
    rhs(end) = 0;
%     rhs(1) = exactSoln(xc(1),t,v);
%     rhs(end) = exactSoln(xc(end),t,v);
    u = A \ rhs;
    
%     Exact = exactSoln(xc,t,v);
    
    figure(1)
    plot(xc, u,'r--.')%,xc,Exact,'b')
    axis([-1,1,-1.5,1.5])
%     pause(0.01)
    drawnow
end

%%
%===Chebyshev Tau===%

%---Parameters---%
v = 0.01;

%---Grid---%
N = 64;
[~,xc] = chebdifmat(N,2,1);
dt = 0.0937;
T = 4;

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

B = I;
B(end-1,end-1) = 0; B(end,end) = 0;

%---Initialization---%
% IC = sin(2*pi*xc) + sin(pi*xc)/2;
IC = sin(pi*xc);
u = DCT*IC;
t = 0;

plot(xc,IC)
axis([-1,1,-1.5,1.5])
pause(0.1)

%---Algorithm---%
while t<T
    t = t+dt;
    u = A \ ( B * ( u - dt*(DCT*( (iDCT*u).*(iDCT*Dspec*u) )) ) );
%     u = A \ ( C * u );
    
    figure(2)
    plot(xc, iDCT*u)
    axis([-1,1,-1.5,1.5])
    drawnow
%     pause(0.01)
end

%%
%===Fourier Galerkin===%

%---Parameters---%
v = 1e-2;

%---Grid---%
N = 64;
dx = 2/N;
x = -1:dx:1-dx; x = x'; %N = N+1;
% x = 0:dx:2*pi-dx; x = x';
T = 10;
dt = 0.01; % did not converge for dt = 1e-4

%---Matrices---%
if mod(N,2)==0
    k = [0:ceil(N/2)-1 0 -floor(N/2)+1:-1];
else
    k = [0:ceil(N/2)-1 -floor(N/2):-1];
end
k = 2*pi*k;

I = eye(N);
A = I + dt*v*diag(k.^2);
B = 1i*diag(k);

%---Initialization---%
IC = sin(pi*x);
u = fft(IC);
t = 0;

%---Algorithm---%
while t<T
    t = t+dt;
    u = A \ ( u - dt*fft( (ifft(u).*ifft(B*u)) ) );
    
    u_phys = ifft(u);
    u_phys = [u_phys ; u_phys(1)];
    
    figure(3)
    plot([x;1],u_phys)
    axis([-1,1,-1.5,1.5])
    drawnow
end