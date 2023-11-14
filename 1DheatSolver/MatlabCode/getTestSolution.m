function testSol = getTestSolution(a,b,D,kx,type,eps,BC)

kt = 3*pi; 
if(strcmp(type,'trig'))
    testSol.ue  = @(x,t)     cos(kx*x)*cos(kt*t);
    testSol.ut  = @(x,t) -kt*cos(kx*x)*sin(kt*t);
    testSol.ux  = @(x,t) -kx*sin(kx*x)*cos(kt*t);
    testSol.uxx = @(x,t) -kx*kx*cos(kx*x)*cos(kt*t);
elseif(strcmp(type,'poly'))
    testSol.ue  = @(x,t) (x.^2 + x + 1)*(t^2 + t + 1);
    testSol.ut  = @(x,t) (x.^2 + x + 1)*(2*t + 1);
    testSol.ux  = @(x,t) (2*x + 1)*(t^2 + t + 1);
    testSol.uxx = @(x,t) 2*(t^2 + t + 1);
elseif(strcmp(type,'exact'))
    % if(strcmp(BC,'D'))
        A = 200;
        B = 100;
    % elseif(strcmp(BC,'N'))
    %     A = 0;
    %     B = 200;
    % end
    C = -1.123;
    testSol.ue  = @(x,t) (A*cos(kx*x + C) + B*sin(kx*x + C))*exp(-D*(kx*kx)*t);
    testSol.ux  = @(x,t) (-kx*A*sin(kx*x + C) + kx*B*cos(kx*x + C))*exp(-D*(kx*kx)*t);
    testSol.ut  = @(x,t) -D*(kx*kx)*testSol.ue(x,t);
    testSol.uxx = @(x,t) -kx*kx*testSol.ue(x,t);
elseif(strcmp(type,'exact1'))
    A = eps/2;
    B = 20;
    C = 10;
    testSol.ue  = @(x,t) A*x.*x + B*x + C + ((2*A)/eps)*t;
    testSol.ux  = @(x,t) 2*A*x + B;
    testSol.ut  = @(x,t) (2*A)/eps;
    testSol.uxx = @(x,t) 2*A;
end

if(strcmp(BC,'D'))
    testSol.ga  = @(t) testSol.ue(a,t);
    testSol.gb  = @(t) testSol.ue(b,t);
elseif(strcmp(BC,'N'))
    testSol.ga  = @(t) -testSol.ux(a,t);
    testSol.gb  = @(t)  testSol.ux(b,t);
end

testSol.u0  = @(x) testSol.ue(x,0);
testSol.f   = @(x,t) testSol.ut(x,t) - D*testSol.uxx(x,t);

end
