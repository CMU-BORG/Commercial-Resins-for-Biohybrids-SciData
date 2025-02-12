function [x,y] = HingeDisplacement(E,x0,y0,theta0,params,Fy,s_samp)

A = params.A;
I = params.I;
L = params.L;
w = params.w;

dt = params.dt;
ds = params.ds;
dds = params.dds;

Mc = params.Mc;
Nc = params.Nc;
n = params.n;

s = linspace(0,1,n)';
rho = 1;
g = 9.8;

% setting up beam equation
DELTA_t = params.DELTA_t;
t = 0;
T = 2;
t_ = t:DELTA_t:T;
T = t_(end); % making sure we get to the maximum force value
i=0;

kappa_0 = diff(theta0)./diff(s); kappa_0 = [kappa_0(1);kappa_0];
kappa_0 = interp1(s,kappa_0,linspace(0,1,n),"linear","extrap");
ksi = 0*x0; eta = 0*y0; eps = 0*x0;

res = zeros(6*n,1);
res_0 = res;

X = res;

% setting initial angle data
i=3;
X(n*(i-1)+(1:n)) = theta0;
X0 = X;

% Solving the dynamics

kk = 0;

while t<=T
    kk = kk + 1;
    fx = 0;
    fy = 0;

    FX_point = 0;
    FY_point = Fy * t / T;

    % Setting up the boundary condition indicators

    non_BCs = 1:6*n;
%     BCs
    BC = [];
    BC_val = @(X) [];
    % no horizontal displacement at left end
    i = 1; %index of the state vector
    alpha = 1; %position index
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X);0*X(n*(i-1)+alpha,:)];
    
    % no vertical displacement at left end
    i = 2; %index of the state vector
    alpha = 1; %position index
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X);0*X(n*(i-1)+alpha,:)];
    
    % no angle at the left end
    i = 3; %index of the state vector
    alpha = 1; %position index
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X); theta0(alpha) + 0*X(n*(i-1) + alpha,:)];
    
    %no shear stress at the right end
    i = 4;
    alpha = n;
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X); FY_point*cos(X(n*(3-1)+alpha,:)) * Nc];%(1000)*(X(n*(2-1) + alpha,:)) ]; %
    
    %no axial stress at the right end
    
    i = 5;
    alpha = n;
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X);FY_point*sin(X(n*(3-1)+alpha,:)) * Nc];
    
    %no moment at the right end
    i = 6;
    alpha = n;
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X); FY_point*cos(X(n*(3-1)+alpha,:))*(w/2) * Mc];
    
    non_BCs(BC) = []; % removing components associated with boundary conditions



    F = @(X) QuasistaticForceBalance(X,fx,fy,0,[dt,ds,E,A,I,L,rho,g,Nc,Mc],dds,x0/ds,y0/ds,theta0);

    F_BC = @(X) [X(BC,:) - BC_val(X);ReturnSubset(F,X,non_BCs)];
    
    F_BC(X0);

    % tol = 1e-6;
    % err = 1e10;
    % Xi = X0;
    % iter = 0;
    % while (err > tol)&(iter<100)
    %     Ji = Jacobian(F_BC,Xi,0.001);
    %     Fi = F_BC(Xi);
    %     dX = Ji \ -Fi;
    % 
    %     if any(isnan(dX))
    %         disp("NaNs")
    %         pause
    %     end
    % 
    %     Xi = Xi + dX;
    %     err = sum(dX.^2);
    %     iter = iter + 1;
    % end

    J_inv = inv(Jacobian(F_BC,X0,1e-4));
    F0 = F_BC(X0);
    
    iter = 0;
    iter_max = 1000;
    tol = 1e-6;
    err = 1e8;
    
    while (err > tol)&(iter<iter_max)
        iter = iter + 1;
    
        % calculate the step size
        dX = - 0.5 * J_inv * F0;
    
        % update solution and calculate new function value
        Xi = X0 + dX;
        Fi = F_BC(Xi);
    
        df = Fi - F0;
    
        err = norm(dX) / n;
    
        % update the inverse Jacobian estimate
        J_inv = J_inv + (1 / ( dX'*J_inv*df ))*( dX - J_inv*df ) * dX' * J_inv;
    
        % re initializing
        X0 = Xi;
        F0 = Fi;
    
    end
    % fprintf("Residual norm: %f \n",norm(Fi))

    X = Xi;
    t = t + DELTA_t;
    X0 = X;
    
    ksi = ds*X(0*n+1:1*n);
    eta = ds*X(1*n+1:2*n);
    dx = dds*(x0/ds + ksi/ds);
    dy = dds*(y0/ds + eta/ds);
    eps = sqrt( ( dx ).^2 + ( dy ).^2 ) - 1; 
    theta = X(2*n+1:3*n);

end

x_ = x0 + ksi;
y_ = y0 + eta;
s_ = linspace(0,1,length(x_));
x = interp1(s_,x_,s_samp);
y = interp1(s_,y_,s_samp);

end