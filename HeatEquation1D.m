% giai phuong trinh truyen nhiet 1 chieu
%       u_t = c^2 u_xx + f(x, t)
%       u(0, t) = p(t)  , t>=0
%       u(L, t) = q(t)  , t>=0
%       u(x, 0) = phi(x), 0<= x <= L

%input  f = f(x, y, t)
%       phi = u(x,0) dieu kien ban dau
%       p = u(0,t), q = u(L,t) dieu kien bien
%       
%Luoi sai phan: 
%       [0, L] khoang khong gian
%       I so luong khoang chia khong gian
%       [0, T] khoang thoi gian
%       N so luong khoang chia thoi gian
%ouput  v: Nghiem so
%       u_ex: nghiem chinh xac

%-------------------------
% Input
f = @(x, t) (1+4*pi*pi)*exp(t)*sin(2*pi*x);
phi = @(x) sin(2*pi*x);
p = @(t) 0;
q = @(t) 0;
L = 1;
I = 20;
T = 1;
N = 800;
c = 1;
u_ex = @(x, t) exp(t)*sin(2*pi*x);
%-------------------------------------
% Phuong phap sai phan huu han
dt = T/N;
dx = L/I;

x = zeros(I+1,1);
t = zeros(N+1,1);

% tinh x_i theo khong gian
for i=1:I
    x(i+1) = i*dx;
end

% tinh t_n theo thoi gian
for n=1:N
    t(n+1) = n*dt;
end

% tinh v
v = zeros(N+1, I+1);
%-------------------------------------
% Thiet lap dieu kien ban dau
v(1, :) = feval(phi, x);

% Thiet lap dieu kien bien
v(:, 1) = feval(p, t);
v(:, I+1) = feval(q, t);

% Dieu kien on dinh
lamda = c^2*dt/dx^2;
%-------------------------------------
% Giai phuong trinh

% Phuong phap hien
v1 = FDM_Explicit_1D(I, N, f, dx, dt, x, t, v, lamda);
% Phuong phap an
v2 = FDM_Implicit_1D(I, N, f, dx, dt, x, t, v, lamda, p, q);
% Phuong phap Crank Nicolson
v3 = CrankNicolson_1D(I, N, f, dx, dt, x, t, v, lamda);
%-------------------------------------
% Tinh nghiem chinh xac
uex = zeros(N+1, I+1);
for n=1:N+1
    for i=1:I+1
        uex(n, i) = u_ex(x(i), t(n));
    end
end
%-------------------------------------
% Tinh nghiem tai thoi diem t co dinh (t = 0.95)
t0 = 0.95;
uex1 = zeros(1, I+1);
v_1 = zeros(1, I+1);
for i=1:I+1
    uex1(i) = u_ex(x(i), t0);
    v_1(i) = v1(N*t0+1, i);
    v_2(i) = v2(N*t0+1, i);
    v_3(i) = v3(N*t0+1, i);
end
%-------------------------------------
%Ve do thi nghiem so

%-------------------------------------
figure(1)
mesh(x, t, v1)
title('Numerical Solution by Explicit FDM')
xlabel('x')
ylabel('t')
zlabel('u')

figure(2)
mesh(x, t, v2)
title('Numerical Solution by Implicit FDM')
xlabel('x')
ylabel('t')
zlabel('u')

figure(3)
mesh(x, t, v3)
title('Numerical Solution by Crank Nicolson')
xlabel('x')
ylabel('t')
zlabel('u')

figure(4)
plot(x, v_1, 'mo', x,v_2,'b',x, v_3, 'g--', x,uex1,'r*-');
title('Exact and numeric solution at time t = 0.95')
xlabel('x'); 
ylabel('u');
legend('Explicit FDM', 'Implicit FDM', 'Crank Nicolson', 'Exact solution');

%----------------------------------
% Phuong phap sai phan hien
function sol1 = FDM_Explicit_1D(I, N, f, dx, dt, x, t, v, lamda)
    % Dieu kien on dinh cua phuong phap
    disp('Luoc do hien');
    if lamda > 1/2
         fprintf("Bai toan co nghiem khong on dinh. ");     
        %return;
    end
    for n = 1:N
       for i = 2:I
           v(n+1, i) = (1-2*lamda)*v(n, i) + ...
                        lamda*(v(n, i-1)+v(n, i+1)) + ...
                        dt*feval(f, x(i), t(n));
       end
    end
    
    sol1 = v;
end
%------------------------------------------
% Phuong phap sai phan an
function sol2 = FDM_Implicit_1D(I, N, f, dx, dt, x, t, v, lamda, p, q)
    disp('Luoc do an');
    a = -lamda;
    d = 1 + 2*lamda;
    % Thiet lap ma tran A là ma tran 3 duong cheo
    A = diag(d*ones(1,I+1)) + diag(a*ones(1,I),1) + diag(a*ones(1,I),-1);
    A(1,1) = 1;
    A(I+1, I+1) = 1;
    A(1, 2) = 0;
    A(I+1, I) = 0;

    % Tinh F
    F = zeros(N+1, I+1);
    F(:, 1) = feval(p, t);
    F(:, I+1) = feval(q, t);    
    for n = 2:N+1
        for i = 2:I
            F(n, i) = v(n-1, i) + dt * feval(f,x(i),t(n));
        end
        v1 = inv(A) * (F(n, :))';
        v(n, :) = v1';
    end
    sol2 = v;
end
%-------------------------------------------
% Phuong phap Crank Nicolson
function sol3 = CrankNicolson_1D(I, N, f, dx, dt, x, t, v, lamda)
    disp('Luoc do Crank Nicolson');
    a = -lamda;
    d1 = 2 + 2*lamda;
    d2 = 2 - 2*lamda;
    A = diag(d1*ones(1,I+1)) + diag(a*ones(1,I),1) + diag(a*ones(1,I),-1);
    A(1,1) = 1;
    A(I+1, I+1) = 1;
    A(1, 2) = 0;
    A(I+1, I) = 0;
    B = diag(d2*ones(1,I+1)) + diag(-a*ones(1,I),1) + diag(-a*ones(1,I),-1);
    B(1,1) = 1;
    B(I+1, I+1) = 1;
    B(1, 2) = 0;
    B(I+1, I) = 0;

    % Tinh F   
    F = zeros(N+1, I+1);
    for n = 1:N
        for i = 2:I
            F(n, i) = 2 * dt * feval(f,x(i),t(n));
        end
        v1 = inv(A) * B * (v(n,:))' + inv(A) * (F(n, :))';
        v(n+1, :) = v1';
    end
    sol3 = v;
end




















