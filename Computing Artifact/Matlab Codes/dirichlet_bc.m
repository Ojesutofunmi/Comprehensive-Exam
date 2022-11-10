%---------------------------------------------------------------------%
%This function computes the dirichlet boundary condition.
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%

function q = dirichlet_bc(qe,q,npoin)
    
    
    q(1) =  qe(1);
    q(npoin) =  qe(npoin);

end