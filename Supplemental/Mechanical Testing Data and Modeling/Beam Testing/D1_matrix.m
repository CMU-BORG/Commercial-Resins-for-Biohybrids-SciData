function [ddx,ddy] = D1_matrix(Nx,Ny,Lx,Ly)
%% returns the derivative matrix for the first order derivatives in both directions
%{
    This function utilizes a central difference method for the central grid
    points and forward/backwards difference method for the boundary grid
    points.
%}

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);

Dx = zeros(Nx,Nx); % d/dx block for one row
% boundary conditions 
Dx(1,1:4) = (1/(6*dx))*[-11, 18, -9, 2];
Dx(2,1:4) = (1/(6*dx))*[-2, -3, 6, -1];
Dx(end-1,end-3:end) = (1/(6*dx))*[1, -6, 3, 2];
Dx(end,end-3:end) = (1/(6*dx))*[-2, 9, -18, 11];
for i = 3:Nx-2
    Dx(i,i-2:i+2) = (1/(12*dx))*[1, -8, 0, 8, -1];
end
[indx_i,indx_j] = find(Dx~=0);
indx_d1 = find(Dx~=0);
Dx_vec = Dx(indx_d1);

Nnz = length(indx_i);
indx_i_total = zeros(Nnz*Ny,1);
indx_j_total = zeros(Nnz*Ny,1);
Dx_vec_total = zeros(Nnz*Ny,1);

for i=1:Ny
    indx_i_total(1+Nnz*(i-1):Nnz*i) = indx_i + Nx*(i-1);
    indx_j_total(1+Nnz*(i-1):Nnz*i) = indx_j + Nx*(i-1);
    Dx_vec_total(1+Nnz*(i-1):Nnz*i) = Dx_vec;
end

ddx = sparse(indx_i_total,indx_j_total,Dx_vec_total,Nx*Ny,Nx*Ny);


indy_i_total = [];
indy_j_total = [];
Dy_vec_total = [];

% inner lower diag
indy_i_total = [indy_i_total,1+Nx:(Ny-1)*Nx];
indy_j_total = [indy_j_total,1:(Ny-2)*Nx];
Dy_vec_total = [Dy_vec_total,-1*ones(size(1+Nx:(Ny-1)*Nx))];
% inner upper diag
indy_i_total = [indy_i_total,1+Nx:(Ny-1)*Nx];
indy_j_total = [indy_j_total,1+2*Nx:Ny*Nx];
Dy_vec_total = [Dy_vec_total,1*ones(size(1+Nx:(Ny-1)*Nx))];

for i=1:Nx
    % starting row
    indy_i_total = [indy_i_total,i*ones(1,3)];
    indy_j_total = [indy_j_total,(1:Nx:3*Nx - 1) + (i-1)];
    Dy_vec_total = [Dy_vec_total, [-3,4,-1]];

    % ending row
    indy_i_total = [indy_i_total,Nx*(Ny-1) + i*ones(1,3)];
    indy_j_total = [indy_j_total,(Nx*(Ny-1):-Nx:Nx*Ny-3*Nx) + i];
    Dy_vec_total = [Dy_vec_total, [3,-4,1]];
end

Dy_vec_total = Dy_vec_total / (2*dy);

ddy = sparse(indy_i_total,indy_j_total,Dy_vec_total,Nx*Ny,Nx*Ny);
