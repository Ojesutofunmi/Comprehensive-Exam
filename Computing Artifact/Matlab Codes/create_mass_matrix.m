
%This function computes the LGL grid and elements.
function element_mass_matrix = create_mass_matrix(intma,coord,Ne,ngl,nq,wnq,L)

%Initialize
element_mass_matrix=zeros(ngl,ngl,Ne);
x=zeros(ngl,1);

% loop over elements
for e=1:Ne
   
   %Store Coordinates
   for i=1:ngl
       I=intma(i,e);
       x(i)=coord(I);
   end
   
   dx=x(ngl)-x(1);
      
   %Do LGL Integration
   for k=1:nq %loop over integration points

      for j=1:ngl %loop over columns
          psi_j=L(j,k);

          for i=1:ngl %loop over rows
              psi_i=L(i,k);
              element_mass_matrix(i,j,e)=element_mass_matrix(i,j,e) + (wnq(k)*(dx/2))*psi_i*psi_j;
         end %i

      end %j

   end %k
   
end %e



      