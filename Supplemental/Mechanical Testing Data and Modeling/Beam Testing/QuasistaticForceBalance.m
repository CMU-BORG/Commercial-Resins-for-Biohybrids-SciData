function res = QuasistaticForceBalance(X,fx,fy,m,params,dds,x0,y0,theta0)
    % returns the force balance residuals which are equal to 0 at the current
    % solution of the backwards euler integration
    
    %% Pulling out model parameters
    n = length(X)/6;

    dt = params(1);
    ds = params(2);
    E = params(3);
    A = params(4);
    I = params(5);
    L = params(6);
    rho = params(7);
    g = params(8);

    Nc = params(9);
    Mc = params(10);

    %% Separating out the current state variables
    ksi = X(0*n+1:1*n,:);
    eta = X(1*n+1:2*n,:);
    theta = X(2*n+1:3*n,:);

    Q = X(3*n+1:4*n,:);
    N = X(4*n+1:5*n,:);
    M = X(5*n+1:6*n,:);

    dx = dds*(x0 + ksi);
    dy = dds*(y0 + eta);
    eps = sqrt( ( dx ).^2 + ( dy ).^2 ) - 1; 

    %% Setting up the residuals equations

    % tangent and normal vector
    n_hat = [dx ./ (eps + 1); dy ./ (eps + 1)];
    t_hat = [-dy ./ (eps + 1); dx ./ (eps + 1)];
    
    % Force and Moment Balance
    r1 = dds*(N.*n_hat(1:n,:)) + dds*(Q.*t_hat(1:n,:)) + (dt^2 / (rho*A*ds) )*fx;
    r2 = dds*(N.*n_hat(n+1:end,:)) + dds*(Q.*t_hat(n+1:end,:)) + (dt^2 / (rho*A*ds) )*fy;
    r3 = dds*M + (A*ds^2/I)*Q ; 
    
    % atan argument
    theta_calc = 0*theta;
    i1 = dx <= 0;
    i2 = dx > 0;
    theta_calc(i1) = mod(atan2( dy(i1) , dx(i1) ),2*pi);
    theta_calc(i2) = atan2( dy(i2) , dx(i2) );

    % Angle constraint and constitutive equations
    r4 = theta - theta_calc; 
    r5 = N - (E*(1 + m)*dt^2 / (rho*ds^2) )*eps ;
    r6 = M - (E*(1 + m)*dt^2 / (rho*ds^2) )*(dds*(theta-theta0) );

    %% Compiling the residual vector
    res = [r1;r2;r3;r4;r5;r6]; 

end

