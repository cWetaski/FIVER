%   AUTHOR: Sebastian Sas-Brunser, Charles Wetaski
%   LAST CHECKED: 2024-06-07 (Charles Wetaski)

function n = getNormalVector(P,VS_logical,ns,Vxyz)
    % GETNORMALVECTOR Returns a normal vector of the specified opaque surface voxel
    % The normal vector in the discrete space is calculated using the algorithm described in [1].
    % [1]: Normal Computation for Discrete Surfaces in 3D Space. Grit Thürmer and Charles A Wüthrich  
    %             (Eurographics 1997, volume 16, number 3)
    %
    % Inputs:
    %   P (1x3 double (int)):       Voxel coordinate of the surface voxel
    %   VS_opaq (3D logical):       Voxel space of opaque voxels
    %   ns (double (int)):          Neighbourhood size (see [1])
    %   Vxyz (1x3 double):          Relative size of voxel in each dimension.
    % Outputs:
    %   n (1x3 double):             Surface normal vector (normalized)

    w=1; % weighting factor (=1 , no weighting)
    
    [sz_x,sz_y,sz_z] = size(VS_logical);
    size_VS = [sz_x,sz_y,sz_z];

    if length(P) == 1 %% assume encoded as linear indices
        [Px,Py,Pz] = ind2sub(size_VS,P);
        P = [Px,Py,Pz]; % replace P with subscripts
    end
    
    if all(P-ns>=1) && all(P+ns<=size_VS) % voxel has enough neighboors
    
        N_p=VS_logical(P(1)-ns:P(1)+ns,P(2)-ns:P(2)+ns,P(3)-ns:P(3)+ns); %neighborhood

    else
        adjacent_boundaries = zeros(2,3); % (1x3 double): row 1 is lower bound, row 2 is upper bound
        a1=P(1)-ns;
        b1=P(1)+ns;
        
        a2=P(2)-ns;
        b2=P(2)+ns;
        
        a3=P(3)-ns;
        b3=P(3)+ns;
        
        if a1<1
            adjacent_boundaries(1,1) = 1-a1;
            a1=1;
        end
        if b1>size_VS(1)
            adjacent_boundaries(2,1) = b1-size_VS(1);
            b1=size_VS(1);
        end
        
        if a2<1
            adjacent_boundaries(1,2) = 1-a2;
            a2=1;
        end
        if b2>size_VS(2)
            adjacent_boundaries(2,2) = b2-size_VS(2);
            b2=size_VS(2);
        end
        
        if a3<1
            adjacent_boundaries(1,3) = 1-a3;
            a3=1;
        end
        if b3>size_VS(3) 
            adjacent_boundaries(2,3) = b3-size_VS(3);
            b3=size_VS(3);
        end
        
        N_p=VS_logical(a1:b1,a2:b2,a3:b3); % neighborhood excluding boundaries
        % Pad array based on boundaries, replicating
        pad_dir = ["pre","post"];
        for ii = 1:2
            N_p = padarray(N_p,adjacent_boundaries(ii,:),"replicate",pad_dir(ii));
        end
    end
    


    [ii,jj,kk]=ind2sub(size(N_p),find(N_p~=0));

    ii=(ii-(ns+1))*Vxyz(1);
    jj=(jj-(ns+1))*Vxyz(2);
    kk=(kk-(ns+1))*Vxyz(3); % converted -ns:ns indexes, scaled by Vxyz

    v=[ii jj kk];

    ii(~any(v,2))=[];
    jj(~any(v,2))=[];
    kk(~any(v,2))=[]; % delete the (i,j,k)=(0,0,0) row.

    v=v(any(v,2),:); % remove the (i,j,k)=(0,0,0) row.

    ss=sqrt(sum(v.^2,2)); % produces column with the sqrt(i^2+j^2+k^2)


    n=[-sum(ii.*w./ss) -sum(jj.*w./ss) -sum(kk.*w./ss)]; % not normalized normal vector at P_hit
    
    n(size_VS == 1) = 0; % if 2D voxel space

    n=n/norm(n); % normalize vector
end %function