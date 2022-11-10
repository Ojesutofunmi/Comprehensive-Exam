

tic

%Input Data
%---------------------------------------------------------------------%

time_final=1; %final time

N = 4;    %Interpolation Order
Ne = 20;  %Number of Elements

integration_points = 1; %=1 for LGL and =2 for LG - LGL is default, but feel free to experiment
integration_type = 2;   %=1 is inexact and =2 is exact - used to compute how many quadrature points we need
space_method_type = 'cg'; %CG or DG - assumes that the code can switch between CG and DG

cfl=0.25; %dt controlled by Courant_max


iplot_solution=1; %Switch to Plot or Not.
iplot_conv=1; 

mu = 0.014;  % diffusivity
u = 2;       % velocity

[q,qe,coord] = Solver_function(N,Ne,u,mu,cfl,time_final,integration_points,space_method_type,integration_type);

%Compute Norm

l2_norm = norm(q-qe,2)/norm(qe,2);
disp(['L2 = ',num2str(l2_norm),' N = ',num2str(N),' Ne = ',num2str(Ne),' npoin = ',num2str(npoin)])

% %Plot Solution
if (iplot_solution == 1)
    figure(1);
    plot_handle=plot(coord,q,'r-');
    set(plot_handle,'LineWidth',2);
    hold on
    plot_handle=plot(coord,qe,'b--');
    set(plot_handle,'LineWidth',2);

    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);

%     title_text=[main_text ': Ne = ' num2str(Ne) ', N = ' num2str(N)', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
%     title([title_text],'FontSize',18);
%     set(gca, 'FontSize', 18);
end


