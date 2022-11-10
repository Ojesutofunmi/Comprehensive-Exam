%---------------------------------------------------------------------%
%This function computes the initial solution.
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%
function [qe] = initial_condition(coord,npoin)

%Initialize
qe=zeros(npoin,1);

%Generate Grid Points
for i=1:npoin
  x=coord(i);  
  qe(i)=exp( -64.0*x^2 ); 
  
end %i      


      
