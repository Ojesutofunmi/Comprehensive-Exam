%Input Data
%---------------------------------------------------------------------%

time_final = 1; %final time

order = [1,2,3,4,8];         % Interpolation Order
NE = [2,4,8,16,32];         % Number of Elements
integration_points = 1;    % = 1 for LGL and = 2 for LG
integration_type = 2;      % = 1 is inexact and = 2 is exact - used to compute how many quadrature points we need
space_method_type ='cg';   % CG or DG 

iplot_solution = 1;        % Switch to Plot or Not.
iplot_conv = 1;

mu = 0.014;  % viscosity
u = 1;   % velocity
cf = 0; %disturbance term
Numerical_convergence(order,NE,u,mu,time_final,integration_points,space_method_type,integration_type,...
                        iplot_solution,iplot_conv)
