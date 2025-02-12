function [x,y] = BeamModel(E,x0,y0,theta0,params,Fy,s_samp)
%{
    Geometrically exact Euler Bernoulli cantilevered beam model. For a
    given dead weight Fy, this function returns the deformed position of
    the beam. The rest position of the beam is captured in the x0,y0 and
    theta0 vectors. The parameters of the beam model are captured in the
    params struct, and the Young's modulus is E. The x,y position of the
    deformed beam are resampled at the locations s_samp.
%}


%% Breaking out the model parameters
A = params.A;   % [m^2] cross sectional area
I = params.I;   % [m^4] 2nd moment of area
L = params.L;   % [m] arc length
w = params.w;   % [m] width

dt = params.dt;     % [s] time step (used for normalizing)
ds = params.ds;     % [m] spatial step
dds = params.dds;   % first derivative finite difference matrix

Mc = params.Mc;     % moment scaling
Nc = params.Nc;     % force scaling
n = params.n;       % number of elements in the beam

s = linspace(0,1,n)';   % arc length parameters of the beams
rho = 1;                % [kg/m^3] density (just in normalization, not used)
g = 9.8;                % [m/s^2] gravitational acceleration

%% Setting up beam equation

% "time" vector of iterative method
DELTA_t = params.DELTA_t;   % step size of iterative method
t = 0;  % starting time
T = 2;  % ending time
t_ = t:DELTA_t:T;   % vector of times
T = t_(end); % making sure we get to the maximum force value


% initializing the residuals vector
res = zeros(6*n,1);
X = res;

% setting initial angle data
i=3;
X(n*(i-1)+(1:n)) = theta0;
X0 = X;

%% Solving the model equations

while t<=T
    fx = 0;
    fy = 0;

    FY_point = Fy * t / T; % scaling the dead weight based on iteration number

    %% Setting up the boundary condition indicators
    non_BCs = 1:6*n;
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
    
    % shear stress at the right end
    i = 4;
    alpha = n;
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X); FY_point*cos(X(n*(3-1)+alpha,:)) * Nc];%(1000)*(X(n*(2-1) + alpha,:)) ]; %
    
    % axial stress at the right end
    
    i = 5;
    alpha = n;
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X);FY_point*sin(X(n*(3-1)+alpha,:)) * Nc];
    
    % moment at the right end
    i = 6;
    alpha = n;
    BC = [BC,n*(i-1)+alpha];
    BC_val = @(X) [BC_val(X); FY_point*cos(X(n*(3-1)+alpha,:))*(w/2) * Mc];
    
    non_BCs(BC) = []; % removing components associated with boundary conditions

    %% Setting up the system of equations
    F = @(X) QuasistaticForceBalance(X,fx,fy,0,[dt,ds,E,A,I,L,rho,g,Nc,Mc],dds,x0/ds,y0/ds,theta0); % force balance equations
    F_BC = @(X) [X(BC,:) - BC_val(X);ReturnSubset(F,X,non_BCs)]; % force balance equations with boundary conditions enforced
       
    %% Solving the system of equations with Broyden's Method
    F0 = F_BC(X0);
    J_inv = inv(Jacobian(F_BC,X0,1e-4)); % initializing the inverse Jacobian
    
    % setting of the solver
    iter = 0;
    iter_max = 1000;
    tol = 1e-6;
    err = 1e8;
    
    % iteration loop
    while (err > tol)&(iter<iter_max)
        iter = iter + 1;
    
        % calculate the step size
        dX = - 0.5 * J_inv * F0;
    
        % update solution and calculate new function value
        Xi = X0 + dX;
        Fi = F_BC(Xi);
        
        % difference in force
        df = Fi - F0;
    
        err = norm(dX) / n;
    
        % update the inverse Jacobian estimate
        J_inv = J_inv + (1 / ( dX'*J_inv*df ))*( dX - J_inv*df ) * dX' * J_inv;
    
        % re initializing
        X0 = Xi;
        F0 = Fi;
    end
    
    % reinitializing
    X = Xi;
    t = t + DELTA_t;
    X0 = X;

end

% pulling out the final displacements and rescaling
ksi = ds*X(0*n+1:1*n);
eta = ds*X(1*n+1:2*n);

%% Reinterpolating the final position of the beam
x_ = x0 + ksi;
y_ = y0 + eta;
s_ = linspace(0,1,length(x_));
x = interp1(s_,x_,s_samp);
y = interp1(s_,y_,s_samp);

end