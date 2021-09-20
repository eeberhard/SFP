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
            %   or figure ID where the graphical editing will be rendered
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
