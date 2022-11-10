%---------------------------------------------------------------------%
%This function computes the exact solution.
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%

function qe = exact_solution(coord,npoin,time,u,mu,Cf)


    %Initialize
    qe=zeros(npoin,1);
    
    %Generate Grid Points
    
    for i=1:npoin
        x = coord(i);
        
        if(mu == 0)
            w = 8;
            xbar = u*time;
            qe(i) = exp(-w^2*(x-xbar)^2);

        elseif(Cf ~= 0)
            qe(i) = sin(2*pi*x); %sin(2*pi);
        else
            w = 2*pi;
            qe(i) = sin(w*x)*exp(u*x/(2*mu))*exp(-(u^2/(4*mu) + w^2*mu)*time);
        
        end
    
    end     
end

      
