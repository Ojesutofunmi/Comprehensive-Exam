%---------------------------------------------------------------------%
%This code contains the function for numerical convergence and plots.
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%
function Numerical_convergence(order,NE,u,nu,time_final,Cf,integration_points,space_method_type,integration_type,...
                        iplot_solution,iplot_conv)

    l2_norm = zeros(length(order),length(NE));
    
    for iN = 1:length(order)
    
        N = order(iN);
    
         cfl = 1/(N+1);
    
        for ie = 1:length(NE)

            
    
            Ne = NE(ie);

            %cfl = 1/(Ne+1);
            [q,qe,coord] = Solver_function(N,Ne,u,nu,cfl,time_final,Cf,integration_points,space_method_type,integration_type);
    
            % Compute Norm
            if(Cf == 0)
                l2_norm(iN,ie) = norm(q-qe,2)/norm(qe,2);
            end
    
        end

        if(Cf == 0)
            
            l2_norm(iN,:);

            fprintf("\n\n");
        end
    
        
    
    end
    
    % Plot Solution
    if (iplot_solution == 1)
        %figure(1);
        plot_handle=plot(coord,q,'r-');
        set(plot_handle,'LineWidth',2);
        hold on
        if(Cf == 0)
            plot_handle=plot(coord,qe,'b--');
        end
        set(plot_handle,'LineWidth',2);
    
        xlabel('x','FontSize',18);
        ylabel('q(x,t)','FontSize',18);
        legend({' CG','Exact'},'Location','northeast')
        title_text=[space_method_type ' method,   integration-type = ' num2str(integration_type)];
        title([title_text],'FontSize',14);
    
    end
    
    % convergence plot
    if (Cf == 0 && iplot_conv == 1)
        figure(2);
        for iN = 1:length(order)
             Np = NE*order(iN)+1; % to calculate the number of points
            if(order(iN)>8)
                p = polyfit(log(Np(1:2)), log(l2_norm(iN,1:2)), 1);
            else
                p = polyfit(log(Np), log(l2_norm(iN,:)), 1);
            end
            txt = ['N = ',num2str(order(iN)), ' rate = ', num2str(round(p(1)))];
            plot_handle = loglog(Np, l2_norm(iN,:), '.-','MarkerSize',18,'DisplayName',txt); hold on
            set(plot_handle,'LineWidth',2);
            xlabel('N_p','FontSize',14);
            ylabel('L^2 Norm','FontSize',14);
            title('CG exact')
            %xticks(NE)
            %xlim([NE(1), NE(end)])
            
            
        end
    
        hold off
        legend show
    
    end

end