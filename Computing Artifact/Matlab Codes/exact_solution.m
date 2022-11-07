
function qe = exact_solution(coord,npoin,time,u,mu)


%Initialize
qe=zeros(npoin,1);
w = pi;

%Generate Grid Points

for i=1:npoin
  x = coord(i);

  qe(i) = cos(w*(x-u*time))*exp(-mu*w^2*time);
  
end     


      
