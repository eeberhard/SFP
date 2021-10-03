classdef SFPGUI < SFP
    %SFPGUI
    %   Graphical editor for the spherical frame projection class, which
    %   allows the reachable regions for each frame axis to be manually
    %   adjusted to account for sparse or noisy samples.
    %
    %   2021 Enrico Eberhard
    %
    %   See also: SFP
    
    properties (Access = private, Hidden)
        h   % Private data container handle
    end
    
    methods
        function obj = SFPGUI(data, resolution)
            % SFPGUI Constructor
            %   Construct a Spherical Frame Projection with a graphical
            %   editor for adjusting the boundary regions of each axis.
            %
            %   This constructor takes the same arguments as the base SFP
            %   class.
            %
            %   See also: SFP
            arguments
                data double
                resolution (1,1) int16 ...
                    {mustBeGreaterThan(resolution, 0), ...
                    mustBeLessThan(resolution, 5)} = 3
            end
            obj = obj@SFP(data, resolution);
            
            msg = ['The SFP GUI allows you to adjust the projected ' ...
                'spherical frame regions by manually selecting or ' ...
                'deselecting vertices that are associated with the ' ...
                'interior area. The X, Y and Z frame regions will be ' ...
                'presented for editing in turn. Close each figure ' ...
                'window to proceed. Finally, the composite SFP will ' ...
                'be plotted.'];
            uiwait(msgbox(msg, 'SFP GUI', 'help'));
            
            for ax = ['X', 'Y', 'Z']
                [~, fig] = obj.edit(ax, figure());
                uiwait(fig)
            end
            
            obj.updatePublicMembers();
            obj.plot();
        end
        
        function [obj, fig] = edit(obj, ax, fig)
            % Edit the SFP graphically
            %
            %   sfp.edit() will open the graphical editor for each axis
            %   in a new figure. Reachable regions can be defined manually
            %   by selecting or deselecting icosphere vertices.
            %   Once the range of motion for a given axis has been edited
            %   to satisfaction, the figure window can be closed.
            %
            %   sfp.edit(ax) will allow editing a single axis region,
            %   where ax is one of 'x', 'y' or 'z'
            %
            %   sfp.edit(..., fig) takes an additional figure handle
            %   or figure ID where the graphical editing will be rendered.
            arguments
                obj
                ax (1, 1) char {mustBeAx(ax, 1)}
                fig (1,1) {mustBeFigure(fig)} = figure()
            end
            if isnumeric(fig)
                fig = round(fig);
            end
            figure(fig);
            clf(fig);
            
            hSphere = patch('Faces', obj.Sphere.Faces, ...
                'Vertices', 0.99*obj.Sphere.Vertices);
            hold on;
            if ~isempty(obj.Tri.(ax))
                hFaces = trisurf(obj.Tri.(ax));
            else
                hFaces = struct();
            end
            
            hFrame = scatter3(obj.Pts.(ax)(:,1), ...
                obj.Pts.(ax)(:,2), obj.Pts.(ax)(:,3), 'ro');
            
            hIdx = scatter3(obj.Sphere.Vertices(obj.Idx.(ax), 1), ...
                obj.Sphere.Vertices(obj.Idx.(ax), 2), ...
                obj.Sphere.Vertices(obj.Idx.(ax), 3), 'k.');
            
            hBoundary = plot3(nan, nan, nan);
            hold off;
            
            title(sprintf('Editing %s Frame Region', upper(ax)));
            
            hSphere.FaceColor = [1, 1, 1];
            hSphere.FaceAlpha = 0.8;
            hSphere.EdgeAlpha = 0.2;
            
            hFaces.EdgeAlpha = 0.3;
            hFaces.FaceAlpha = 0.5;
            hFaces.FaceColor = [0.5 0.5 0.5];
            hFaces.HitTest = 'off';
            hFaces.PickableParts = 'none';
            
            hFrame.MarkerFaceColor = 'flat';
            hFrame.MarkerFaceAlpha = 0.25;
            hFrame.PickableParts = 'none';
            
            hIdx.SizeData = 100;
            hIdx.HitTest = 'off';
            hIdx.PickableParts = 'none';
            
            hBoundary.Color = 0.25 * hFaces.FaceColor;
            hBoundary.LineWidth = 2.5;
            hBoundary.HitTest = 'off';
            hBoundary.PickableParts = 'none';
            
            axis equal
            axis vis3d
            
            switch upper(ax)
                case 'X'
                    view(90, 0)
                case 'Y'
                    view(180, 0)
                case 'Z'
                    view(180, 90)
            end
            
            obj.h.(ax) = struct(...
                'hSphere', hSphere, ...
                'hFaces', hFaces, ...
                'hFrame', hFrame, ...
                'hIdx', hIdx, ...
                'hBoundary', hBoundary);
            
            hSphere.ButtonDownFcn = @(~, hit)(obj.onClick(hit, ax));
            obj.redrawBoundary(ax);
        end
        
        function obj = rotate(obj, R)
            % Rotate the sampled data and regions
            %
            %   sfp.rotate(R) will rotate the stored orientation samples,
            %   the underlying icosphere and the generated surface regions
            %   by a supplied rotation R, which can take the form of either
            %   a homogenous 3x3 rotation matrix or a unit quaternion.
            %
            %   Example:
            %       % create an SFP from some orientation set Q
            %       % and edit the regions with the GUI
            %       sfp = SFPGUI(Q);
            %       % rotate the SFP by 60 degrees about the X axis
            %       sfp.rotate([0.866, 0.5, 0, 0]);
            %       % plot the newly rotated SFP
            %       sfp.plot()
            arguments
                obj
                R double {validateRotation(R)} = eye(3)
            end
            
            [q, R] = validateRotation(R);
            
            for t = 1:size(obj.q, 1)
                obj.q(t, :) = SFP.hamilton(q, obj.q(t, :)) .* [1 -1 -1 -1];
                obj.Pts.X(t, :) = obj.rotateVec(obj.q(t,:), [1 0 0]);
                obj.Pts.Y(t, :) = obj.rotateVec(obj.q(t,:), [0 1 0]);
                obj.Pts.Z(t, :) = obj.rotateVec(obj.q(t,:), [0 0 1]);
            end
            obj.Sphere.Vertices = obj.Sphere.Vertices * R';
            obj.triangulateRegions();
            obj.updatePublicMembers();
        end
    end
    
    methods (Access = protected, Hidden)
        function obj = onClick(obj, hit, ax)
            I = obj.findNearest(obj.Sphere.Vertices, hit.IntersectionPoint);
            
            f = find(obj.Idx.(ax) == I, 1);
            if ~isempty(f)
                obj.Idx.(ax)(f) = [];
            else
                obj.Idx.(ax)(end+1) = I;
            end
            
            obj.h.(ax).hIdx.XData = obj.Sphere.Vertices(obj.Idx.(ax), 1);
            obj.h.(ax).hIdx.YData = obj.Sphere.Vertices(obj.Idx.(ax), 2);
            obj.h.(ax).hIdx.ZData = obj.Sphere.Vertices(obj.Idx.(ax), 3);
            
            obj.triangulateRegions();
            if ~isempty(obj.Tri.(ax))
                obj.h.(ax).hFaces.Faces = obj.Tri.(ax).ConnectivityList;
                obj.h.(ax).hFaces.Vertices = obj.Tri.(ax).Points;
            end
            
            obj.redrawBoundary(ax);
        end
        
        function obj = redrawBoundary(obj, ax)
            data = nan(1, 3);
            for b = obj.Boundaries.(ax)
                data = [data; obj.Sphere.Vertices(b{1}, :); nan(1, 3)]; %#ok<AGROW>
            end
            
            obj.h.(ax).hBoundary.XData = data(:, 1);
            obj.h.(ax).hBoundary.YData = data(:, 2);
            obj.h.(ax).hBoundary.ZData = data(:, 3);
        end
    end
    
end

%% Local functions

function q = quatFromRot(R)
%%QUATFROMROT
%   Quaternion from rotation matrix
%
%   q = quatFromRot(R) takes a 3x3 rotation matrix and returns a quaternion
%       for that rotation.
%
%   Code adapted from:
%   http://www.peterkovesi.com/matlabfns/Rotations/matrix2quaternion.m
%
%   2018 Enrico Eberhard
arguments
    R (3,3) double
end
if any(isnan(R))
    q = nan(1,4);
    return
end

% Find rotation axis as the eigenvector having unit eigenvalue
% Solve (R-I)v = 0;
[v,d] = eig(R-eye(3));

% The following code assumes the eigenvalues returned are not necessarily
% sorted by size. This may be overcautious.
d = diag(abs(d));   % Extract eigenvalues
[~, ind] = sort(d); % Find index of smallest one
if d(ind(1)) > 0.001   % Hopefully it is close to 0
    warning('Rotation matrix is dubious');
end

axis = v(:,ind(1)); % Extract appropriate eigenvector
if abs(norm(axis) - 1) > .0001     % Debug
    warning('Non-unit rotation axis');
end

% Now determine the rotation angle
twocostheta = trace(R)-1;
twosinthetav = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
twosintheta = axis'*twosinthetav;

theta = atan2(twosintheta, twocostheta);

q = [cos(theta/2); axis*sin(theta/2)]';

end

function R = rotFromQuat(q)
%%ROTFROMQUAT
%   Rotation matrix from quaternion
%
%   q = rotFromQuat(R) takes a unit quaternion and returns a 3x3 rotation
%   matrix for that rotation.
%
%   2018 Enrico Eberhard
arguments
    q (1,4) double
end
if any(isnan(q))
    R = nan(3, 3);
    return
end

w = q(1);
x = q(2);
y = q(3);
z = q(4);
R = [1 - 2*y*y - 2*z*z, 2*x*y - 2*w*z, 2*x*z + 2*w*y;
    2*x*y + 2*w*z, 1 - 2*x*x - 2*z*z, 2*y*z - 2*w*x;
    2*x*z - 2*w*y, 2*y*z + 2*w*x, 1 - 2*x*x - 2*y*y];

end


%% Custom input validation functions
function mustBeFigure(fig)
if ~(isa(fig, 'matlab.ui.Figure') || (isnumeric(fig) && fig >= 1))
    msg = sprintf('Input type \"%s\" cannot be converted to a figure.', class(fig));
    throwAsCaller(MException('SFP:mustBeFigure', msg));
end
end

function mustBeAx(ax, dim)
arguments
    ax
    dim = 0
end
if dim == 1
    if sum(upper(ax) == 'XYZ') ~= 1
        msg = sprintf("Ax input must be 'X', 'Y', or 'Z'.");
        throwAsCaller(MException('SFP:mustBeAx', msg));
    end
else
    if ~all(sum(upper(ax) == ('XYZ')'))
        msg = sprintf("Ax input must only contain 'X', 'Y', or 'Z'.");
        throwAsCaller(MException('SFP:mustBeAx', msg));
    end
end
end

function [q, R] = validateRotation(R)
if size(R, 1) == 3 && size(R, 2) == 3
    % rotation matrix
    q = quatFromRot(R);
elseif size(R, 1) == 1 && size(R, 2) == 4
    q = R;
    R = rotFromQuat(q);
elseif size(R, 1) == 4 && size(R, 2) == 1
    q = R';
    R = rotFromQuat(q);
else
    throwAsCaller(MException('SFP:mustBeRotation', ...
        'Input must be a quaternion or 3x3 rotation matrix'));
end
end