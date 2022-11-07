
function Numerical_convergence(order,NE,u,mu,time_final,integration_points,space_method_type,integration_type,...
                        iplot_solution,iplot_conv)

    l2_norm = zeros(length(order),length(NE));
    
    for iN = 1:length(order)
    
        N = order(iN);
    
         cfl = 1/(N+1);
    
        for ie = 1:length(NE)

            Ne = NE(ie);

            %cfl = 1/(Ne+1);
            [q,qe,coord] = Solver_function(N,Ne,u,mu,cfl,time_final,integration_points,space_method_type,integration_type);
    
            % Compute Norm
    
            l2_norm(iN,ie) = norm(q-qe,2)/norm(qe,2);
    
        end
        l2_norm(iN,:)
    
        
    
    end
    
    % Plot Solution
    if (iplot_solution == 1)
        figure(1);
        plot_handle=plot(coord,q,'r-');
        set(plot_handle,'LineWidth',2);
        hold on
        plot_handle=plot(coord,qe,'b--');
        set(plot_handle,'LineWidth',2);
    
        xlabel('x','FontSize',18);
        ylabel('q(x,t)','FontSize',18);
    
    end
    
    % convergence plot
    if (iplot_conv == 1)
        figure(2);
        for iN = 1:5
            p = polyfit(log(NE), log(l2_norm(iN,:)), 1);
            txt = ['N = ',num2str(order(iN)), ' rate = ', num2str(p(1))];
            plot_handle = loglog(NE, l2_norm(iN,:), '.-','MarkerSize',18,'DisplayName',txt); hold on
            set(plot_handle,'LineWidth',2);
            %xticks(NE)
            %xlim([NE(1), NE(end)])
            
            
        end
    
        hold off
        legend show
    
    end

end