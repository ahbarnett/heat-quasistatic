function [un,t] = quasiStatic(i1,i2,ia,ib,Ngx,Nt,dx,dt,ga,gb,eps)

I = i1:i2;

% prepare the matrix for the first elliptic solve
A= sparse(Ngx, Ngx);

for i=I
    A(i, i-1)= 1/(dx*dx);
    A(i, i)  = -2.0/(dx*dx);
    A(i, i+1)= 1/(dx*dx);
end

A(ia, ia)= 1;
A(ib, ib)= 1;

w0 = zeros(Ngx, 3); 
w1 = zeros(Ngx, 3); 

cnt = 1; 
for n = Nt-1:1:Nt+1
    t = n*dt;

    rhs = zeros(1,Ngx)';
    rhs(ia) = ga(t);
    rhs(ib) = gb(t);
    w0(:,cnt) = A\rhs;

    rhs = w0(:,cnt); 
    rhs(ia) = 0;
    rhs(ib) = 0;
    w1(:,cnt) = A\rhs; % what is the initial condition here? 

    cnt = cnt + 1; 
end

v0 = w0(:,2); 
v1 = (1/(2*dt))*(w1(:,3) - w1(:,1)); 
un = v0 + eps*v1; 

t = Nt*dt; 
end
