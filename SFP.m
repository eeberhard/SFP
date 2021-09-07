classdef SFP < handle
    %SFP Spherical Frame Projection
    %   Use the Spherical Frame Projection method to generate and
    %   visualize a the range of motion for a 3D rotational joint.
    %
    %   2021 Enrico Eberhard
    %
    %   See also: SFPGUI
    
    properties (SetAccess = protected)
        N (1, 1) int16  % Icosphere order (resolution)
        q (:, 4) double % Underlying orientation data (array of quaternions)
        Sphere (1,1) struct % Icosphere face and vertex structure
        X  % Structure of faces, vertices and boundaries for the X region
        Y  % Structure of faces, vertices and boundaries for the Y region
        Z  % Structure of faces, vertices and boundaries for the Z region
    end
    properties (Access = protected, Hidden)
        Boundaries % IDs of the Sph.Vertices
        Faces      % IDs of the Sph.Faces that are counted as "interior"
        Idx        % IDs of the Sph.Vertices that match the
        Pts        % Frame points
        Tri        % triangulation
    end
    
    methods
        function obj = SFP(data, resolution)
            %SFP Constructor
            %   obj = SFP(data) will construct an sfp object from
            %   unit quaternion orientation data in the form of an
            %   Nx4 matrix, where N is the number of samples and each
            %   row contains the quaternion coefficients as w, x, y, z.
            %
            %   obj = SFP(data, resolution) will accept an integer
            %   resolution argument between 1 and 5 (default 3) to specify
            %   the order (number of faces) on the icosphere
            arguments
                data double
                resolution (1,1) int16 ...
                    {mustBeGreaterThan(resolution, 0), ...
                    mustBeLessThan(resolution, 5)} = 3
            end
            if size(data, 2) == 4
                obj.q = data;
            elseif size(data, 1) == 4
                obj.q = data';
            else
                error('SFP:invalidArgument', ...
                    'Input data must be size 4 in one dimension');
            end
            
            obj.N = resolution;
            obj.init();
            
            obj.indexFramePoints();
            obj.triangulateRegions();
            obj.updatePublicMembers();
        end
        
        function obj = fillHoles(obj, threshold, axes)
            % Fill holes in projected regions
            %   sfp.fillHoles() will try to identify interior holes
            %   in the projected region. This is useful when the sampled
            %   orientation data has insufficent coverage, particularly at
            %   higher resolutions of the SFP icosphere. It returns a
            %   modified SFP object with the identified holes filled in by
            %   one level.
            %
            %   fillHoles(threshold) takes an integer threshold value
            %   which refers to the number of edges an interior hole must
            %   have before it is filled. A larger threshold will result in
            %   more aggressive infilling.
            %
            %   fillHoles(..., axes) takes a character array of axes for
            %   the regions to fill (for example, 'X' or 'YZ'). The default
            %   behaviour is to fill all regions ('XYZ').
            %
            %   The call to fillHoles may need to be repeated, as the
            %   infilling adds one layer of faces to the enclosing regions.
            %   If the radius of the hole therefore is more than two faces
            %   at any point, a smaller hole may remain.
            %
            %   Note that the projected region for a given frame axis
            %   may contain "true" holes.
            %
            %   See also: findHoles, trimIslands
            arguments
                obj
                threshold (1,1) int16 = 20
                axes (1, :) char {mustBeAx(axes)} = 'XYZ'
            end
            for ax = upper(axes)
                holes = obj.findHoles(ax, threshold);
                
                if isempty(holes)
                    continue
                end
                
                connected = [];
                for hole = holes
                    for vertex = hole{1}'
                        connected = [connected; ...
                            obj.findConnectedVertices(obj.Sphere.Faces, vertex)]; %#ok<AGROW>
                    end
                end
                
                obj.Idx.(ax) = unique([obj.Idx.(ax); connected]);
            end
            
            obj.triangulateRegions();
        end
        
        function obj = trimIslands(obj, threshold, axes)
            % Trim islands in projected regions
            %   sfp.trimIslands() will try to identify exterior
            %   islands in the projected region. This is useful when the
            %   sampled orientation data is noisy and has some outliers,
            %   particularly at higher resolutions of the SFP icosphere.
            %   It returns a modified SFP object with the identified
            %   islands trimmed by one level.
            %
            %   trimIslands(threshold) takes an integer threshold value
            %   which refers to the number of edges an island must
            %   have before it is trimmed. A larger threshold will result in
            %   more aggressive trimming.
            %
            %   trimIslands(..., axes) takes a character array of axes for
            %   the regions to trim (for example, 'X' or 'YZ'). The default
            %   behaviour is to trim all regions ('XYZ').
            %
            %   The call to trimIslands may need to be repeated, as the
            %   trimming removes one layer of faces from the exterior regions.
            %   If the radius of the island therefore is more than two faces
            %   at any point, a smaller island may remain.
            %
            %   Note that the projected region for a given frame axis must
            %   be continguous to represent range of motion of a spherical
            %   joint with no external factors. Any islands indicate poor
            %   fitting of the data.
            %
            %   See also: findHoles, fillHoles
            arguments
                obj
                threshold (1,1) int16 = 20
                axes (1, :) char {mustBeAx(axes)} = 'XYZ'
            end
            
            for ax = upper(axes)
                [~, islands] = obj.findHoles(ax, threshold);
                
                if isempty(islands)
                    continue
                end
                
                for island = islands
                    for vertex = island{1}'
                        obj.Idx.(ax)(obj.Idx.(ax) == vertex) = [];
                    end
                end
            end
            
            obj.triangulateRegions();
        end
        
        function [holes, islands] = findHoles(obj, ax, threshold)
            % Find holes and islands in a projected region
            %   holes = sfp.findHoles(ax) will return a cell array of
            %   vertex indices corresponding to identified holes in the
            %   specific region. The input ax can be one of 'X', 'Y' or 'Z'.
            %
            %   holes = sfp.findHoles(ax, threshold) takes an integer
            %   threshold value to determine the hole size cut-off.
            %
            %   [holes, islands] = sfp.findHoles(...) will also return
            %   a cell array of vertex indices identifying exterior island
            %   regions, using the same axis and threshold for both.
            %
            %   See also: fillHoles, trimIslands
            arguments
                obj
                ax (1, 1) char {mustBeAx(ax, 1)}
                threshold (1,1) int16 = 20
            end
            
            [holes, islands] = deal({});
            ax = upper(ax);
            for b_ = obj.Boundaries.(ax)
                b = b_{1};
                if numel(b) >= threshold
                    continue
                end
                
                % march along the boundary and accumulate the ratio of
                % connected faces to boundary faces
                nConnected = 0;
                nEmpty = 0;
                for id = 1:numel(b) - 1
                    for face = obj.findConnectedFaces(obj.Sphere.Faces, b(1))'
                        if ismember(obj.Sphere.Faces(face, :), obj.Faces.(ax), 'rows')
                            nConnected = nConnected + 1;
                        else
                            nEmpty = nEmpty + 1;
                        end
                    end
                end
                
                % if each boundary vertex is generally connected to more
                % existing faces than empty faces, they surround a hole
                if nConnected > nEmpty
                    holes{end+1} = b; %#ok<AGROW>
                else
                    islands{end+1} = b; %#ok<AGROW>
                end
            end
        end
        
        function plot(obj, fig, axes)
            % PLOT Plot the SFP
            %   sfp.plot() renders the spherical frame project in a new
            %   figure. It draws the faces of the projected regions
            %   corresponding to the reachability of the X, Y, Z frame axes
            %   in red, green and blue, respectively.
            %   The border of each region is also drawn.
            %
            %   plot(fig) takes an integer or existing figure handle on
            %   which to draw the SFP. The 'hold' property of the figure is
            %   set to 'on', and the SFP is drawn with a unit radius,
            %   centered on the origin, with XYZ frame axes aligned.
            %
            %   plot(..., axes) takes a character array of axes for
            %   the regions to draw (for example, 'X' or 'YZ'). The default
            %   behaviour is to draw all regions ('XYZ').
            arguments
                obj
                fig (1,1) {mustBeFigure(fig)} = figure()
                axes (1,:) char {mustBeAx(axes)} = 'XYZ'
            end
            
            if isnumeric(fig)
                fig = round(fig);
            end
            figure(fig);
            
            % associate a color with each region
            cols = struct('X', [1 0.5 0.5], ...
                'Y', [0.5 1 0.5], ...
                'Z', [0.5 0.5 1]);
            
            % offset the height of each region slightly
            % to properly render overlapping areas
            scale = struct('X', 0.99, ...
                'Y', 1.00, ...
                'Z', 1.01);
            
            hold on;
            hSphere = patch('Faces', obj.Sphere.Faces, ...
                'Vertices', 0.98*obj.Sphere.Vertices);
            
            hSphere.FaceColor = [1, 1, 1];
            hSphere.FaceAlpha = 0.2;
            hSphere.EdgeAlpha = 0.05;
            
            for ax = upper(axes)
                tri = obj.Tri.(ax);
                tri = triangulation(tri.ConnectivityList, ...
                    scale.(ax) * tri.Points);
                h = trisurf(tri);
                
                h.EdgeAlpha = 0.3;
                h.FaceAlpha = 0.5;
                h.FaceColor = cols.(ax);
                
                
                hFrame = scatter3(obj.Pts.(ax)(:,1), ...
                    obj.Pts.(ax)(:,2), obj.Pts.(ax)(:,3), 'k.');
                
                for b = obj.Boundaries.(ax)
                    p = plot3(scale.(ax) * obj.Sphere.Vertices(b{1}, 1), ...
                        scale.(ax) * obj.Sphere.Vertices(b{1}, 2), ...
                        scale.(ax) * obj.Sphere.Vertices(b{1}, 3));
                    p.Color = 0.75 * cols.(ax);
                    p.LineWidth = 2.5;
                end
                
            end
            hold off;
            axis equal
            axis vis3d
            axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);
        end
    end
    
    methods (Access = protected, Hidden)
        
        function obj = init(obj)
            % Initialise object properties
            obj.Sphere = icosphere(obj.N);
            
            for ax = ['X', 'Y', 'Z']
                obj.Boundaries.(ax) = {};
                obj.Faces.(ax) = [];
                obj.Pts.(ax) = [];
            end
        end
        
        function obj = indexFramePoints(obj)
            % Find the closest vertices to the orientation frame vectors
            %   ...
            [I, obj.Pts.X, obj.Pts.Y, obj.Pts.Z] = ...
                deal(zeros(size(obj.q, 1), 3));
            
            for t = 1:size(obj.q, 1)
                obj.Pts.X(t, :) = obj.rotateVec(obj.q(t,:), [1 0 0]);
                obj.Pts.Y(t, :) = obj.rotateVec(obj.q(t,:), [0 1 0]);
                obj.Pts.Z(t, :) = obj.rotateVec(obj.q(t,:), [0 0 1]);
                
                % find closest V to p
                I(t, 1) = obj.findNearest(obj.Sphere.Vertices, obj.Pts.X(t, :));
                I(t, 2) = obj.findNearest(obj.Sphere.Vertices, obj.Pts.Y(t, :));
                I(t, 3) = obj.findNearest(obj.Sphere.Vertices, obj.Pts.Z(t, :));
            end
            obj.Idx.X = unique(I(:, 1));
            obj.Idx.Y = unique(I(:, 2));
            obj.Idx.Z = unique(I(:, 3));
        end
        
        function obj = triangulateRegions(obj)
            % Fit faces and boundaries to the vertices of each region
            %   ...
            warnstate = warning;
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
            for ax = ['X', 'Y', 'Z']
                obj.Faces.(ax) = obj.findInteriorFaces(obj.Sphere.Faces, ...
                    obj.Idx.(ax));
                
                if isempty(obj.Faces.(ax))
                    obj.Tri.(ax) = [];
                    obj.Boundaries.(ax) = {};
                    continue
                end
                obj.Tri.(ax) = triangulation(obj.Faces.(ax), ...
                    obj.Sphere.Vertices);
                
                I = obj.Tri.(ax).freeBoundary();
                
                gaps = [0, find(I(2:end, 1) ~= I(1:end-1, 2))', size(I, 1)];
                
                boundaries = cell(1, numel(gaps) - 1);
                
                for g = 1:(numel(gaps) - 1)
                    boundaries{g} = [I((gaps(g) + 1):gaps(g + 1), 1); ...
                        I(gaps(g + 1), 2)];
                end
                
                obj.Boundaries.(ax) = boundaries;
            end
            warning(warnstate);
        end
        
        function obj = updatePublicMembers(obj)
            % Update the public data structures for each frame axis
            %   ...
            for ax = ['X', 'Y', 'Z']
                obj.(ax) = struct('Faces', obj.Faces.(ax), ...
                    'Vertices', obj.Sphere.Vertices, ...
                    'Boundaries', obj.Boundaries.(ax), ...
                    'FramePoints', obj.Pts.(ax));
            end
        end
    end
    
    methods (Static, Access = protected)
        function I = findNearest(A, B)
            [~, I] = min(vecnorm(A - B, 2, 2));
        end
        
        function faces = findConnectedFaces(F, vertex)
            % for the vertex ID, return all connected face IDs
            faces = unique(mod(find(F == vertex), size(F, 1)));
        end
        
        function vertices = findConnectedVertices(F, vertex)
            % for the vertex ID, return all connected vertex IDs
            faces = SFP.findConnectedFaces(F, vertex);
            v = F(faces, :);
            vertices = unique(v(:));
        end
        
        function FI = findInteriorFaces(F, I)
            % for a list of faces and vertex IDs, return all faces that are
            % exclusively connected to the given vertices
            FI = zeros(0, 3);
            for f = F'
                % check that all corners are inside
                if any(f(1) == I) && any(f(2) == I) && any(f(3) == I)
                    FI = [FI; f']; %#ok<AGROW>
                end
            end
        end
        
        function v_ = rotateVec(q, v)
            % rotate a vector by a unit quaternion
            arguments
                q (1,4) double
                v (1,3) double
            end
            q = q ./ norm(q);
            
            V = hamilton(q, hamilton([0 v], q.*[1 -1 -1 -1]));
            
            v_ = V(2:4);
            
            function q = hamilton(q2, q1)
                q(1) = q1(1)*q2(1)-dot(q1(2:4),q2(2:4));
                q(2:4) = q1(1)*q2(2:4) + q2(1)*q1(2:4) ...
                    + cross(q1(2:4),q2(2:4));
            end
        end
        
    end
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

