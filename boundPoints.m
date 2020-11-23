function boundPoints(q, N)

if N < 1 || N > 4
    error('N must be between 1 and 4');
end

% max edge length of unit icosphere with degree N
EDGE_LENGTH = sqrt(2) * [0.5465, 0.2759, 0.1383, 0.0692, 0.0346];

FV = icosphere(N);
F = FV.Faces;
V = FV.Vertices;
% VN = FV.VertexNormals;

% for each orientation, find closest vertex to axis vector
[I, Iy, Iz] = deal(zeros(1, length(q)));
for t = 1:length(q)
    X = quatMultVec(quatConj(q(t,:)), [1 0 0]);
    Y = quatMultVec(quatConj(q(t,:)), [0 1 0]);
    Z = quatMultVec(quatConj(q(t,:)), [0 0 1]);
    
    % find closest V to p
    I(t) = findNearest(V, X);
    [~, Iy(t)] = min(vecnorm(V - Y, 2, 2));
    [~, Iz(t)] = min(vecnorm(V - Z, 2, 2));
end

I1 = unique(I);
B = findBoundaryPoints(V, I1, EDGE_LENGTH(N));
insideFaces = findInteriorFaces(F, I1);

% march around boundary

% ve = V(E,:);
% vne = VN(E,:);
% B = ve(1, :);
% while true
%     next = findNearest(ve, ve(1,:));
%     B = [B; ve(next, :)];
%     
%     D = B(end, :) - B(end - 1, :);
%     
% end


fig = gcf;
clf(fig)

h = patch('Faces', insideFaces, 'Vertices', V);
h.FaceAlpha = 0.5;
h.FaceColor = [1 1 1];
h.EdgeAlpha = 0.2;
axis equal
axis vis3d

hold on;
scatter3(V(I1, 1), V(I1, 2), V(I1, 3), 'LineWidth', 2);
scatter3(V(B, 1), V(B, 2), V(B, 3), 'LineWidth', 5);
hold off;

end

%% Subfunctions

function I = findNearest(A, B)
    [~, I] = min(vecnorm(A - B, 2, 2));
end

function B = findBoundaryPoints(V, I, L)
    B = [];
    for ind = I
        nn = sum(vecnorm(V - V(ind, :), 2, 2) < L);
        n = sum(vecnorm(V(I, :) - V(ind, :), 2, 2) < L);
        if n < nn && n > 2
            B = [B ind]; %#ok<AGROW>
        end
    end
end

function FI = findInteriorFaces(F, I)
    FI = zeros(0, 3);
    for f = F'
        % check that all corners are inside
        if any(f(1) == I) && any(f(2) == I) && any(f(3) == I)
            FI = [FI; f']; %#ok<AGROW>
        end
    end
end

function findBoundaryEdges %#ok<DEFNU>
% GOAL: order the boundary points by connected
% edges so that a boundary line can be drawn 

% 0: create an index list containing the first boundary point

% 1: get the last boundary point from the index list

% 2: find all (interior) faces that contain this point

% 3: for each face, find all other boundary points connected
% to the original point by an edge

% 4: discard any boundary points that are already in the index list

% 5: add the remaining boundary point to the index list.
% a. If there is exactly one candidate, add it to the index list.
% b. If there are multiple candidates, choose the boundary point
% belonging to a face containing only boundary points.
% c. If there are no candidates, but still unindexed boundary points,
% end the index list, start a new list with the first next
% unindexed boundary point

% 7: if there are any remaining unindexed boundary points, repeat 
% from one, else end

end