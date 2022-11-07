%This code computes the matrix of the lagrange basis derivative 
% we have L and dL matrix

function [L,dL] = LagrangeMaxtrix(P,Q,xgl, Xq)

    L = zeros(P,Q);
    dL = zeros(P,Q);
    
    
    for k = 1:Q
      
       %Construct Basis
       for i=1:P

          Li = 1;
          dLi = 0;
          
          for j = 1:P
             
             %Basis
             if (j ~= i)
                 
                Li = Li*(Xq(k)-xgl(j))/(xgl(i)-xgl(j));
             
                 prod = 1;

                 if (j ~= i)
                    for l = 1:P
                        
                       %Derivative
                       if (l ~=i && l ~=j)
                          prod = prod*(Xq(k)-xgl(l))/(xgl(i)-xgl(l));
                       end
                       
                    end
                    
                    dLi = dLi + prod/(xgl(i)-xgl(j));
                    
                 end
             end
          end
          
          L(i,k) = Li;
          dL(i,k) = dLi;
       end
    end

end
