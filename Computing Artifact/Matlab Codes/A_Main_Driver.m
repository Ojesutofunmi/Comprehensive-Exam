%---------------------------------------------------------------------%
%This code contains the main driver for the convergence study
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%

time_final = 1; %final time

clc;

%case 1 
%order = [2,3,4,6,8];         % Interpolation Order
%NE = [8,16,32,64];         % Number of Elements     

%case 2 
order = [2,3,4,8,10,12,16];  % Interpolation Order
NE = [4,8,16,32];         % Number of Elements

integration_points = 1;    % = 1 for LGL and = 2 for LG
integration_type = 2;      % = 1 is inexact and = 2 is exact - used to compute how many quadrature points we need
space_method_type ='cg';   % CG or DG 

iplot_solution = 0;        % Switch to Plot or Not.
iplot_conv = 1;

nu = 0.21;  % diffusivity %0.21
u = 2;   % velocity
Cf = 0; %rate of snowfall

Numerical_convergence(order,NE,u,nu,time_final,Cf,integration_points,space_method_type,integration_type,...
                        iplot_solution,iplot_conv)


