function Le = Laplace_matrix(coord,intma,Ne,ngl,nq,wnq,dpsi)
    
    
    Le = zeros(ngl,ngl,Ne);
    x=zeros(ngl,1);

    for e=1:Ne
    
        %Store Coordinates
        for i=1:ngl
           I=intma(i,e);
           x(i)=coord(I);
        end
        
        dx=x(ngl)-x(1);
        
        for k = 1:nq                  % loop through integration points
            wk = wnq(k);                        % quadrature weight
            
            for j = 1:ngl              % loop through columns
                xj = dpsi(j,k);            % value of the derivative of the basis function
                
                for i =1:ngl          % loop through rows
                    xi = dpsi(i,k);        % value of the derivative of the basis function
                    
                    Le(i,j,e) = Le(i,j,e) - (2/dx)*wk*xi*xj;
                end
            end
        end
    end
end
