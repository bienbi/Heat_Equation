% giai pt truyen nhiet hai chieu voi bien Dirichlet 
%input 
%         f = f(x,y,t) 
%         p1 = u(0, y, t)
%         p2 = u(a, y, t)
%         q1 = u(x, 0, t)
%         q2 = u(x, b, t)
%         phi = u(x, y, 0)
%         c he so
%Luoi sai phan
%         [0, a]: khoang chia x
%         [0, b]: khoang chia y
%         T khoang thoi gian
%         I so khoang chia theo x
%         J so khoang chia theo y
%         N so khoang chia theo thoi gian t
%output
%         v nghiem xap xi

%----------------------------------------
% Input
f = @(x,y,t) cos(x+y+t) + 2*sin(x+y+t);
p1 = @(y,t) sin(y+t);
p2 = @(y,t) sin(y+t+1);
q1 = @(x,t) sin(x+t);
q2 = @(x,t) sin(x+t+1);
phi = @(x,y) sin(x+y);
a = 1;
b = 1;
c = 1;
T = 1;
I = 50;
J = 50;
N = 50;
uex = @(x,y,t) sin(x+y+t);

dx = a/I;
dy = b/J;
dt = T/N;

x = zeros(I+1,1);
y = zeros(J+1,1);
t = zeros(N+1,1);
% Tinh x_i
for i=0:I
    x(i+1) = i*dx;
end
% Tinh y_j 
for i=0:J
    y(i+1) = i*dy;
end
% Tinh t_n
for i=1:N
    t(i+1) = i*dt;
end

v = zeros(I+1, J+1, N+1);
% Tinh dieu kien ban dau  
for i=1:I+1
    for j=1:J+1
        v(i, j, 1) = feval(phi, x(i), y(j));
    end
end

% Dieu kien on dinh
lamda1 = c^2*dt/dx^2;
lamda2 = c^2*dt/dy^2;

% Giai phuong trinh

v1 = FDM_Explicit_2D(I, J, N, f, dx, dy, dt, x, y, t, v, lamda1, lamda2, p1, p2, q1, q2);
v2 = FDM_Implicit_2D(I, J, N, f, dx, dy, dt, x, y, t, v, lamda1, lamda2, p1, p2, q1, q2);
v3 = CrankNicolson_2D(I, J, N, f, dx, dy, dt, x, y, t, v, lamda1, lamda2);
%-------------------------------------
% Tinh nghiem chinh xac
u = zeros(I+1, J+1);
for I=1:I+1
    for j=1:J+1
        u(i, j) = uex(x(i), y(j), 0.78);
    end
end
%-------------------------------------
%Ve do thi nghiem so
%-------------------------------------
figure(1)
surf(x, y, v1)
title('Numerical solution by Explicit FDM 3D view')
xlabel('x')
ylabel('y')
zlabel('u')

figure(2)
surf(x, y, v2)
title('Numerical solution by Explicit FDM 3D view')
xlabel('x')
ylabel('y')
zlabel('u')

figure(4)
plot(x, v_1, 'mo', x,v_2,'b',x, v_3, 'g--', x,uex1,'r*-');
title('Exact and numeric solution at time t = 0.95')
xlabel('x'); 
ylabel('u');
legend('Explicit FDM', 'Implicit FDM', 'Crank Nicolson', 'Exact solution');
%----------------------------------
function sol1 = FDM_Explicit_2D(I, J, N, f, dx, dy, dt, x, y, t, v, lamda1, lamda2, p1, p2, q1, q2)
    % Dieu kien on dinh cua phuong phap
    disp('Luoc do hien');
    if lamda1 + lamda2 > 1/2
         fprintf("Bai toan co nghiem khong on dinh. \n");     
        %return;
    end
    % Tinh v
    for n=2:N+1
        for j=1:J+1
            v(1,j,n) = feval(p1,y(j),t(n-1));
            v(I+1,j,n) = feval(p2,y(j),t(n-1));
        end
        for i=1:I+1
            v(i,1,n) = feval(q1,x(i),t(n-1));
            v(i,J+1,n) = feval(q2,x(i),t(n-1));
        end
        for i=2:I
            for j=2:J
                v(i,j,n) = (1-2*lamda1-2*lamda2)*v(i,j,n-1)...
                    +lamda1*(v(i+1,j,n-1)+v(i-1,j,n-1))...
                    +lamda2*(v(i,j-1,n-1)+v(i,j+1,n-1))...
                    +dt*feval(f,x(i),y(j),t(n-1)); 
            end
        end
    end
    sol1 = v(:, :, N+1);
end
%------------------------------------------
function sol2 = FDM_Implicit_2D(I, J, N, f, dx, dy, dt, x, y, t, v, lamda1, lamda2, p1, p2, q1, q2)
    disp('Luoc do an');
    % Thiet lap ma tran A
    A1 = 0.5*lamda1*ones(J-1,1);
    B1 = A1;
    A1(1) = 0;
    B1(I-1) = 0;
    C1 = (1+lamda1)*ones(J-1,1);

    A2 = 0.5*lamda2*ones(I-1,1);
    B2 = A2;
    A2(1) = 0;
    B2(J-1) = 0;
    C2 = (1+lamda2)*ones(I-1,1);


    for n=1:N
        % giai bai toan theo x
        v_temp = zeros(I+1,J+1);
        for j=1:J+1
            v_temp(1, j) = feval(p1, y(j), t(n)+dt/2);
            v_temp(I+1, j) = feval(p2, y(j), t(n)+dt/2);
        end
        for i=1:I+1
            F = zeros(J-1,1);
            for j=2:J
                F(j-1) = 0.5*lamda2*v(i,j-1,n) + (1-lamda2)*v(i,j,n) + 0.5*lamda2*v(i,j+1,n) + 0.5*dt*feval(f,x(i),y(j),t(n)+dt/2);
                if j==2
                    F(j-1) = F(j-1) + v_temp(1,i)*0.5*lamda2;
                end
                if j==J
                    F(j-1) = F(j-1) + v_temp(I+1,i)*0.5*lamda2;
                end
            end
            v_temp(i,:) = GiaiHePT3DuongCheo(J-1,A1,B1,C1,F,v_temp(1,i),v_temp(I+1,i));
        end
        % giai theo y
        v_temp2 = zeros(I+1,J+1);
        for i=1:I+1
            v_temp2(i,1) = feval(q1,x(i),t(n+1));
            v_temp2(i,J+1) = feval(q2,x(i),t(n+1));
        end
        for j=1:J+1
            G = zeros(I-1,1);
            for i=2:I
                G(i-1) = 0.5*lamda1*v(i-1,j) + (1-lamda1)*v(i,j) + 0.5*lamda1*v(i+1,j) + 0.5*dt*feval(f,x(i),y(j),t(n)+dt/2);
                if i==2
                    G(i-1) = G(i-1) + v_temp2(j,1)*0.5*lamda1;
                end
                if i==J
                    G(i-1) = G(i-1) + v_temp2(j,J+1)*0.5*lamda1;
                end
            end
            v_temp2(:, j) = GiaiHePT3DuongCheo(I-1,A2,B2,C2,G,v_temp2(j,1),v_temp2(j,J+1));
        end
        v(:, :, n+1)=v_temp2;
    end
    sol2 = v(:, :, N+1);
end

%-------------------------------------------
function sol3 = CrankNicolson_2D(I, J, N, f, dx, dy, dt, x, y, t, v, lamda1, lamda2)
    disp('Luoc do Crank Nicolson');
    
    sol3 = v(:, :, N+1);
end
%------------------------------------------
function Y  = GiaiHePT3DuongCheo( N,A,B,C,F,y0,yn)
    alpha = zeros(N,1);
    beta = zeros(N,1);
    Y = zeros(N+2,1);
    alpha(1 )= B(1)/C(1);
    beta(1) = F(1)/C(1);
    for i=2:N
        alpha(i)=B(i-1)/( C(i-1)-A(i-1)*alpha(i-1) );
        beta(i)=(A(i-1)*beta(i-1) + F(i-1)) / (C(i-1) - A(i-1)*alpha(i-1));
    end
    %disp(beta);
    %disp(alpha);
    Y(N+1)=(A(N)*beta(N)+F(N))/(C(N)-A(N)*alpha(N));
    for i=1:N-1
        Y(N+1-i) = alpha(N+1-i) * Y(N+2-i) + beta(N+1-i);
        %fprintf('%f %f %f\n',alpha(N+1-i),Y(N+2-i),beta(N+1-i));

    end
    Y(1) = y0;
    Y(N+2) = yn;
end


















