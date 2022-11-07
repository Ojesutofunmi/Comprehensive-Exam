
function [coord,intma,npoin] = create_grid(N,nelem,xgl,space_method_type)

    %Set some constants
    
    ngl = N+1;
    npoin = N*nelem + 1; %number of points for CG 

    xmin=-1;
    xmax=+1;
    dx=(xmax-xmin)/nelem;
    coord=zeros(npoin,1);
    
    %Generate COORD and INTMA for CG
    ip=1;
    coord(1)=xmin;
    for e=1:nelem
       x0=xmin + (e-1)*dx;
       intma(1,e)=ip;
       for i=2:ngl
          ip=ip + 1;
          coord(ip)=(xgl(i)+1 )*dx/2 + x0;
          intma(i,e)=ip;
       end %i
    end %e
    
    if strcmp(space_method_type,'dg')
        
        npoin_dg = ngl*nelem; %number of points for DG
        coord_dg=zeros(npoin_dg,1);

        %Generate COORD and INTMA for DG
        ip=0;
        for e=1:nelem
            for i=1:ngl
                ip=ip+1;
                intma_dg(i,e)=ip;
            end
        end
        
        for e=1:nelem
           for i=1:ngl
              ip_cg=intma(i,e);
              ip_dg=intma_dg(i,e);
              coord_dg(ip_dg)=coord(ip_cg);
           end 
        end

        intma = intma_dg;
        coord = coord_dg;
        npoin = npoin_dg;
    end
end
    