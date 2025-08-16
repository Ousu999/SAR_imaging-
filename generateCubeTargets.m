function [positions, rcs] = generateCubeTargets(center, dims, pointsPerSide, rcsPerPoint)
% generateCubeTargets creates scatter points for a cube or box target
% center       - 3×1 [x; y; z] center position (m)
% dims         - 3×1 [Lx; Ly; Lz] dimensions (m)
% pointsPerSide - number of scatter points per edge
% rcsPerPoint  - RCS value assigned to each scatter point (m²)

    Lx = dims(1); Ly = dims(2); Lz = dims(3);

    % Grid for each face
    linX = linspace(-Lx/2, Lx/2, pointsPerSide);
    linY = linspace(-Ly/2, Ly/2, pointsPerSide);
    linZ = linspace(-Lz/2, Lz/2, pointsPerSide);

    [X, Y] = meshgrid(linX, linY);
    [Xz, Z] = meshgrid(linX, linZ);
    [Yz, Z2] = meshgrid(linY, linZ);

    faces = [];

    % Top & bottom
    faces = [faces; [X(:), Y(:),  Lz/2*ones(numel(X),1)]]; % top
    faces = [faces; [X(:), Y(:), -Lz/2*ones(numel(X),1)]]; % bottom

    % Front & back
    faces = [faces; [X(:),  Ly/2*ones(numel(X),1), Z(:)]]; % front
    faces = [faces; [X(:), -Ly/2*ones(numel(X),1), Z(:)]]; % back

    % Left & right
    faces = [faces; [ Lx/2*ones(numel(Yz),1), Yz(:), Z2(:)]]; % right
    faces = [faces; [-Lx/2*ones(numel(Yz),1), Yz(:), Z2(:)]]; % left

    % Shift to cube center
    positions = (faces' + center);

    % Assign constant RCS to each point
    rcs = rcsPerPoint * ones(1, size(positions, 2));
end
