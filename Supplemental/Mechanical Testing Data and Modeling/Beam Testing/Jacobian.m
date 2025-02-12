function J = Jacobian(F,x,h)
    
    dx = h*ones(size(x));
    x1 = repmat(x,size(x'))+diag(dx);
    x0 = repmat(x,size(x'))-diag(dx);
    dF = F(x1) - F(x0);
    J =  dF ./ (2*dx');

%     tol = 0.01*median(abs(J),"all");
%     J(abs(J)<tol) = 0;
%     J = sparse(J);

end
