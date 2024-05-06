function conn_counts = countConnectivity(A)
    %COUNTCONNECTIVITY Returns 6-connectivity of A
    %   This function counts number of nonzero 6-connected neighbours of each non-ero value of A.
    %   6-connected refers only to voxels which are connected face-to-face. Voxels which have a value of 0 are 
    %   returned as 0, even if they have nonzero neighbours. 
    %   Boundaries are considered nonzero.
    % INPUTS:
    %   A (3D double): Any 3D matrix
    % OUTPUTS:
    %   conn_counts (3D double (int)) [-]: Corresponding 6-connectivity of each nonzero value of A.
    
    A = logical(A); % convert to binary array (i.e. 0 if empty, 1 if not empty)

    [xlim,ylim,zlim] = size(A); % get limits

    B = ones(xlim,ylim,zlim,6); % create 4D array for shifted matrices. Change to from ones to zeros for zero boundaries.

    B(1:(xlim-1),:,:,1) = A(2:xlim,:,:); % shifted -1 along x axis
    B(2:(xlim),:,:,2) = A(1:(xlim-1),:,:); % shifted +1 along x axis
    B(:,1:(ylim-1),:,3) = A(:,2:ylim,:); % shifted -1 along y axis
    B(:,2:ylim,:,4) = A(:,1:(ylim-1),:); % shifted +1 along y axis
    B(:,:,1:(zlim-1),5) = A(:,:,2:zlim); % shifted -1 along z axis
    B(:,:,2:zlim,6) = A(:,:,1:(zlim-1)); % shifted +1 along z axis

    conn_counts = sum(B~=0,4); % sum along 4th dimension to get 3D matrix of original size.
    conn_counts(A == 0) = 0; % set voxels with zero value to zero, since we only want connectivity of nonzero voxels.
end