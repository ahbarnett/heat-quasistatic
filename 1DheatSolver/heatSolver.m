% Solve the heat equation on the unit interval
% Use Rothe's method or the Quasistatic method
% Use a sample exact solution for testing
% This file performs a grid refinement study and records the order of
% accuracy

clear; clf; close all; fontSize = 14;  lineWidth = 2;

%% Parameters
method = 'Quasistatic'; % choose 'Rothe' or 'Quasistatic'
BC     = 'D';           % choose the type of BC 'D' for Dirichlet or 'N' for Neumann

if(~(strcmp(BC,'D')||strcmp(BC,'N')))
    error('BC type does not exist');
end
if(strcmp(BC,'N')&& strcmp(method,'Quasistatic'))
    error('Neumann BC is not supported for Quasistatic method');
end

eps = 0.00001; 
D  = 1/eps;       % the diffusion coefficient
kx = sqrt(eps);   % to use in the exact solution
a  = 0.0;
b  = 1.0;
tFinal = 1.0;
plotOption = 0;

if(strcmp(BC,'N')); BCtype = 'Neumann'; else; BCtype = 'Dirichlet'; end
fprintf('Using the %s method with %s BC and eps = %1.1e\n', method,BCtype,eps);

%% Get sample solution for testing
testSolType = 'exact';
testSol = getTestSolution(a,b,D,kx,testSolType,eps,BC);

numResolutions = 10;
errMax = zeros(1,numResolutions);
%% Begin a grid resoltuion study
for m = 1:numResolutions

    % prepare the grid in space
    Nx = 10*2^(m-1); % number of space intervals

    if(strcmp(BC,'D'))
        ng = 0;
    elseif(strcmp(BC,'N'))
        ng = 1;
    end

    ia = ng + 1;       % index for the grid point at the left boudnary
    ib = ia + Nx;     % index for the grid point at the right boundary
    Ngx= ib + ng;     % total number of grid points including ghost points
    dx = (b-a)/Nx;
    x = zeros(Ngx,1);  % values of x over the grid
    for( ix=1:Ngx )
        x(ix)=a + (ix-ia)*dx;
    end
    i1 = ia+1; % first interior point
    i2 = ib-1; % last interior point

    if(strcmp(BC,'D'))
        I = i1:i2;
    elseif(strcmp(BC,'N'))
        I = ia:ib;
    end

    % allocate space for the solution
    unp1 = zeros(Ngx, 1); % holds U_i^{n+1}

    % time step
    dt = dx; % time step (adjusted below)
    Nt = round(tFinal/dt); % number of time-steps
    dt= tFinal/Nt; % adjust dt to reach tFinal exactly

    % initialize the solution at tn and t(n-1) time levels
    unm1 = testSol.u0(x); % initial conditions
    un   = unm1 + (dt/(eps*dx*dx))*(testSol.u0(x+dx) - 2*unm1 + testSol.u0(x-dx)) + dt*testSol.f(x,0);

    % Apply the time-stepping method
    if(strcmp(method,'Rothe'))
        [un,t] = Rothe_CD2(D, dx, dt, Nt, Ngx, un, unm1,ia,ib,I,testSol.ga,testSol.gb,testSol.f,x,BC);
    elseif(strcmp(method,'Quasistatic'))
        [un,t] = quasiStatic(i1,i2,ia,ib,Ngx,Nt,dx,dt,testSol.ga,testSol.gb,eps);
    elseif(strcmp(method,'BDF'))
       [un,t] = BDF_CD2(D, dx, dt, Nt, Ngx, un, unm1,ia,ib,I,testSol.ga,testSol.gb,testSol.f,x,BC); 
    else
        error('Method name does not exist');
    end

    % Error Analysis
    ue = testSol.ue(x,t);
    err = (un - ue);
    errMax(m) = max(abs(err));

    fprintf('t=%10.4e: Nx=%3d Nt=%4d dt=%9.3e maxErr=%8.2e',t, Nx, Nt, dt, errMax(m));
    if(m==1); fprintf('\n'); else; fprintf(' order=%8.2e\n', log2(errMax(m-1)/errMax(m))); end

    % plot results
    if(plotOption == 1)
        figure(1);
        plot(x,un,'r-o',x,ue,'k-','Linewidth',lineWidth);
        legend('computed', 'true');
        title(sprintf('%s Method t=%1.1e (Nx=%d) dt=%1.1e',method, t, Nx, dt))
        xlabel('x'); ylabel('u'); set(gca, 'FontSize', fontSize); grid on;

        % plot error
        figure(2);
        plot(x,err,'b-x', 'Linewidth', lineWidth);
        legend('Error');
        title(sprintf('%s Method Error t=%1.1e dt=%1.1e',method, t, dt));
        xlabel('x'); ylabel('u'); set(gca, 'FontSize', fontSize); grid on;

        pause;
    end
end