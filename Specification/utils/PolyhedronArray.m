classdef PolyhedronArray < Polyhedron
    %PolyhedronArray Extends the class Polyhedron to contain multiple
    % polyhedrons and to correct the method "contains"

    properties
        PArray = [] % Array of polyhedrons
    end

    methods
        function obj = PolyhedronArray(PArray)
            %PolyhedronArray Constructor
            assert(size(PArray, 2) == 1, "Polyhedron arrays with multiple columns are not supported.")
            obj.PArray = PArray;
        end

        function labels = labelPoints(obj, points)
            %labelPoints Labels points according to their containment in
            % the polyhedrons in the array, e.g., [2, 3, 1] if the points
            % 1, 2, and 3 are contained in the polyhedrons in array ...
            % rows 2, 3, and 1, respectively.
            % Raises an error if labels are ambiguous!

            % Check containment
            inP = contains(obj, points);
        
            % Check if containment unambigous
            temp = sum(inP, 1);
            assert(max(temp) <= 1, "Some points are contained in multiple polyhedrons.")

            % Generate unambigous labelling
            % 1) Uniquely contained points
            [~, labels] = max(inP, [], 1);
            % 2) Points not contained by any polytope
            idx_nc = (temp == 0);
            labels(idx_nc) = 0;
        end

        function inP = contains(obj, points)
            %contains Check whether the polyhedrons contain the points in
            % points. Returns a vector of size nPolyhedrons times nPoints,
            % where an entry in row i and column j indicates that point j
            % is contained in polyhedron i.

            % Switch between 1D and multidimensional case
            dim = size(points, 1);
            nPoints = size(points, 2);
            nPolyhedrons = size(obj.PArray, 1);
            if dim == 1
                %% 1D case
                % Init labels
                inP = ones(nPolyhedrons, nPoints);

                % Iterate over polyhedrons
                for i = 1:nPolyhedrons
                    % Check which points are contained in polyhedron i
                    inP(i, :) = (points > min(obj.PArray(i).V)) & ...
                        (points < max(obj.PArray(i).V));
                end
            elseif dim >= 1
                error("Multi dimensional case: To be written.")
            else
                error("No points given.")
            end
        end
    end
end