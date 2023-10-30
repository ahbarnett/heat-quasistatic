function [un,t] = Rothe_CD2(D, dx, dt, Nt, Ngx, un, unm1,ia,ib,I,ga,gb,f,x, BC)

% prepare the matrix for the elliptic solve
gamma = (2/3)*dt*D;
A= sparse(Ngx, Ngx);

for i=I
    A(i, i-1)= 1/(dx*dx);
    A(i, i)  = -2.0*(1/(dx*dx)) - 1/gamma;
    A(i, i+1)= 1/(dx*dx);
end

if(strcmp(BC,'D'))
    A(ia, ia)= 1;
    A(ib, ib)= 1;
elseif(strcmp(BC,'N'))
    i = ia-1;
    A(i, ia-1)=  1/(2*dx);
    A(i, ia+1)= -1/(2*dx);
    i = ib+1;
    A(i, ib-1)= -1/(2*dx);
    A(i, ib+1)=  1/(2*dx);
end

rhs = zeros(1,Ngx)';

for n = 2:Nt
    t = n*dt;
    % prepare the rhs for the elliptic solve
    for i = I
        rhs(i) = -(4/(3*gamma))*un(i) + (1/(3*gamma))*unm1(i) - (2/(3*gamma))*dt*f(x(i),t);
    end

    if(strcmp(BC,'D'))
        rhs(ia) = ga(t);
        rhs(ib) = gb(t);
    else
        rhs(ia-1) = ga(t);
        rhs(ib+1) = gb(t);
    end

    unp1 = A\rhs;
    unm1 = un;
    un = unp1;
end

end