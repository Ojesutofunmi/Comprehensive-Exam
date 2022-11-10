%---------------------------------------------------------------------%
%This function computes the Differentiation matrix.
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%
function element_differentiation_matrix = create_diff_matrix(ngl,nq,wnq,psi,dpsi)

    %Initialize
    element_differentiation_matrix=zeros(ngl,ngl);
    
    
    %LGL Integration
    for k=1:nq %loop over integration points
      wk=wnq(k);
    
      for j=1:ngl %loop over columns
          psi_j=psi(j,k);
    
          for i=1:ngl %loop over rows
    
              dpsi_i=dpsi(i,k);
              
    
              element_differentiation_matrix(i,j)=element_differentiation_matrix(i,j) + wk*dpsi_i*psi_j;
          end %i
    
      end %j
    
    end %k      
end
