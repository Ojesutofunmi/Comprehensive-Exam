


function [M,D,L] = Matrix_DSS(Me,De,intma,ngl,Ne,npoin,Le)

%Form Global Matrices
M=zeros(npoin); %mass matrix
D=zeros(npoin); %derivative matrix
L = zeros(npoin); % laplacian matrix

for e=1:Ne %loop over elements


    for j=1:ngl %loop over columns of M(e)
        J = intma(j,e); %point to global gridpoint number

        for i=1:ngl %loop over rows of M(e)
            I = intma(i,e); %point to global gridpoint number

            M(I,J)=M(I,J) + Me(i,j,e);
            D(I,J)=D(I,J) + De(i,j);

            L(I,J) = L(I,J) + Le(i,j,e);

        end %i

    end %j


end %e


      