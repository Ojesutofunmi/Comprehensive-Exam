
function [q,qe,coord] = case3(N,Ne,u,mu,cfl,time_final,integration_points,space_method_type,integration_type)
    

    ngl = N+1;
    %Compute Interpolation Points
    [xgl,wgl]=legendre_gauss_lobatto(ngl); %LGL interpolation points
    
    %Compute Integration Points
    if (integration_points == 1) %LGL integration points
        integration_text='LGL';
    if (integration_type == 1) %inexact
        noq = N;
    elseif (integration_type == 2) %exact
        noq=N+1;
    end
    nq=noq + 1;
    [Xq,wnq]=legendre_gauss_lobatto(nq);
    elseif (integration_points == 2) %LG integration points
        integration_text='LG';
        noq = N;
        nq=noq + 1;
        [Xq,wnq]=legendre_gauss(nq);
    end
    
    %Compute Lagrange Polynomial and derivatives
    
    [L,dL] = LagrangeMaxtrix(ngl,nq,xgl,Xq);
    
    %Create Grid - different arrays for CG and DG
    [coord,intma,npoin] = create_grid(N,Ne,xgl,space_method_type);
    
    
    %Choose time-step such that it does not violate cfl but also
    %divides evenly into the final time
    
    if(u ~= 0 && mu == 0)
        dx=coord(2)-coord(1);
        dx = 1/npoin;
        dt=cfl*dx/u;
        ntime=round(time_final/dt);
        dt=time_final/ntime;

    elseif(u == 0 && mu ~= 0)
        dx = coord(2)-coord(1);
        dx = 1/npoin;
        dt=cfl*dx^2;
        ntime=round(time_final/dt);
        dt=time_final/ntime;

    else
        dx = coord(2)-coord(1);
        dx = 1/npoin;
        dt = cfl*min(dx/u , dx^2);
        ntime = round(time_final/dt);
        dt = time_final/ntime;

    end
    
    disp(['dt = ',num2str(dt),' ntime = ',num2str(ntime),' time_final = ',num2str(time_final)])
    
    %main_text=[space_method_type ': ' integration_text];
    
    
    % Create Element Mass, Differentiation and Laplacian matrices
    
    % Create element mass matrix Me
    Me = create_mass_matrix(intma,coord,Ne,ngl,nq,wnq,L);
    
    % Create element differentiation matrix De
    
    De = create_diff_matrix(ngl,nq,wnq,L,dL);
    
    % Create element laplacian matrix Le
    
    Le = Laplace_matrix(coord,intma,Ne,ngl,nq,wnq,dL);
    
    % Form Global Matrices (N_p x N_p)
    
    [Mmatrix,Dmatrix,Lmatrix] = Matrix_DSS(Me,De,intma,ngl,Ne,npoin,Le);
    
    % Create elements boundary fluxes
    if strcmp(space_method_type,'cg') 
        Fmatrix = zeros(npoin,npoin);
    elseif strcmp(space_method_type,'dg')
        if (flux_type == 1)
            Fmatrix = Fmatrix_centered_flux(intma,Ne,npoin,ngl);
        elseif (flux_type == 2)
            Fmatrix = Fmatrix_upwind_flux(intma,Ne,npoin,ngl);
        end
    end
    
    % Create right-hand-size matrix Rmatrix
    % Derive Rmatrix from the matrix form of governing equation 
    
    Rm= u*(Dmatrix - Fmatrix) + mu*Lmatrix;
    
    
    % For BCs dirichlet
%     Mmatrix(1,1) = 1;
%     Mmatrix(npoin,npoin) = 1;

     %Mmatrix(1,:) = 0;
     Mmatrix(1,1) = 1;
     %Mmatrix(npoin,:) = 0;
     Mmatrix(npoin,npoin) = 1;
    
    Rm(1,:) = 0;
    Rm(1,1) = 1;
    Rm(npoin,:) = 0;
    Rm(npoin,npoin) = 1;
    
    %Left-Multiply by Inverse Mass Matrix
    Rmatrix = Mmatrix\Rm;
    
    %Compute Initial Condition
    time=0;
    qe = exact_solution(coord,npoin,time,u,mu);
    
    %Initialize the solution with the initial solution qe
    q=qe;
    
    %Time Integration
    [q,time] = time_integration(q,Rmatrix,time,Ne,dt,coord,u,mu);
    q = (inv(Mmatrix - Lmatrix))*(Mmatrix*F - (Fmatrix*dqdx));
    %Compute Exact Solution
    qe = exact_solution(coord,npoin,time,u,mu);

end