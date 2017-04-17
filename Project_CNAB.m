%% Project - Burgers Equation with CNAB

%%
%===Chebyshev Collocation, CNAB===%

%---Parameters---%
v = 0.00097;

%---Grid---%
N = 64;
[Dm,xc] = chebdifmat(N,2,1);
T = 2;
% dt = 1e-5;
dt = 1e-5;

%---Matrices for CNAB---%
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
    
    plot(xc, u,'ro--',xc,Exact,'b')
    axis([-1,1,0,1.25])
    pause(0.01)
end

%%
%===Chebyshev Tau, CNAB===%

%---Parameters---%
v = 8.7*1e-4;

%---Grid---%
N = 64;
[~,xc] = chebdifmat(N,2,1);
T = 2;
dt = 1e-5;

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
    
    plot(xc, iDCT*u,'ro--',xc,Exact,'b')
    axis([-1,1,0,1.25])
%     pause(0.01)
    drawnow
end

%%
%===Fourier Galerkin, CNAB===%

%---Parameters---%
v = 0.028;

%---Grid---%
N = 64;
dx = 2*pi/N;
x = 0:dx:2*pi-dx; x = x';
T = 4;
dt = 1e-2;

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
    
    plot([x;2*pi], u_phys,'ro--',[x;2*pi], Exact,'b')
    axis([0,2*pi,-3,3])
%     pause(0.01)
    drawnow
end