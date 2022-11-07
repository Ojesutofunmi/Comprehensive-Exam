
function [q0,time] = time_integration(q0,Dhat,time,ntime,dt,coord,u,mu)
    
    
    %Initialize RK coefficients

    RKA = [0, -567301805773 / 1357537059087, -2404267990393 / 2016746695238, ...
            -3550918686646 / 2091501179385, -1275806237668 / 842570457699];
    
    RKB = [1432997174477 / 9575080441755, 5161836677717 / 13612068292357, 1720146321549 / 2090206949498,...
                3134564353537 / 4481467310338, 2277821191437 / 14882151754819];
    
    Npoin = length(q0);
    dq = zeros(Npoin,1);
    qp=q0;
    stages = length(RKA);
    
    %Time Integration
    for itime=1:ntime
        time=time + dt;
        % Compute exact solution
        qe = exact_solution(coord,Npoin,time,u,mu);

        %RK Stages
        for s = 1:stages
            %Create RHS vector
            R=Dhat*qp; %only valid for CG
            %Solve System
            for I=1:Npoin
                dq(I) = RKA(s)*dq(I) + dt*R(I);
                qp(I) = qp(I) + RKB(s)*dq(I);
            end
%             qp = dirichlet_bc(qe,qp,Npoin);
        end %s

        

%         R=Dhat*q0; 
%         qp = q0 +dt*R;

        qp = dirichlet_bc(qe,qp,Npoin);

        %Update Q
        q0 = qp;
    
    end
end
