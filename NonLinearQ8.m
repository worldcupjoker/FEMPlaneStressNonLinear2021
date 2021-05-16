% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: NonLinearQ8.m
% This file solves a 2 demensional non linear plane stress problem using
% finite element method with 8-node quadrilateral element type.
% Input: [file name]

classdef NonLinearQ8 < handle
    properties (Access = private)
        
        % main function run time.
        time;
        
        % Max iteration.
        nMax = 30;
        
        % err
        err = 10 ^ (-7);
        
        % Mode. 0 for tangent, 1 for secant.
        
        % Constants
        E0;
        v0;
        q0;
        
        % nodes - 2D array with 3 columns - node no., x-coordinate, y-coordinate
        % and rows = no. of nodes.
        nodeGlobal;
        totalNodes;
        
        % elements - 2D array with columns - element no., nodes in element and rows
        % = no. of elements .
        elements;
        totalElements;
        
        % Matrices
        K_global;
        f_n;
        f_s;
        u_n;
        f_g;
        relativeResidual = [];
        
        % This is a 4Elements x 6 matrix. First column: x; second column: y;
        % thrid column: sigma_xx; fourth column: sigma_yy; fifth column:
        % sigma_xy; sixth column: sigma_e, von Mises stress.
        % Each element has four flux vectors, with one at each quadrature
        % point.
        stress;
    end
    
    properties (Access = private)
        % Temperary object.
        tempObj;
        
        % Get shape functions for two-point Gaussian quadrature.
        N_a = getShapeFuncQ8(-1 / 3 ^ 0.5, -1 / 3 ^ 0.5, true);
        N_b = getShapeFuncQ8(1 / 3 ^ 0.5, -1 / 3 ^ 0.5, true);
        N_c = getShapeFuncQ8(1 / 3 ^ 0.5, 1 / 3 ^ 0.5, true);
        N_d = getShapeFuncQ8(-1 / 3 ^ 0.5, 1 / 3 ^ 0.5, true);
        
        N_real_a;
        N_real_b;
        N_real_c;
        N_real_d;
        
        
        % Get G matrix (the gradient of shape functions) values for two-point Gaussian Quadrature.
        % The order is a, b, c, d, counter clockwise, starting from the bottom left corner.
        G_a = getShapeFuncQ8(-1 / 3 ^ 0.5, -1 / 3 ^ 0.5, false);
        G_b = getShapeFuncQ8(1 / 3 ^ 0.5, -1 / 3 ^ 0.5, false);
        G_c = getShapeFuncQ8(1 / 3 ^ 0.5, 1 / 3 ^ 0.5, false);
        G_d = getShapeFuncQ8(-1 / 3 ^ 0.5, 1 / 3 ^ 0.5, false);
    end
    
    methods
        
        % Constructor
        function obj = NonLinearQ8(fileName, v_0, E_0, q_0, mode)
            
            % Read in data.
            strGlobalNodes = '*Node';
            strElementNodes = '*Element, type=CPS8R';
            obj.nodeGlobal = readinp(strGlobalNodes, fileName);
            obj.totalNodes = size(obj.nodeGlobal, 1);
            obj.elements = readinp(strElementNodes, fileName);
            obj.totalElements = size(obj.elements, 1);
            obj.v0 = v_0;
            obj.E0 = E_0;
            obj.q0 = q_0;
            
            if (strfind(fileName, 'Prob1') == 1)
                obj.tempObj = Prob1Related(fileName);
            elseif ((strfind(fileName, 'Prob2') == 1))
                obj.tempObj = Prob2Related(fileName);
            elseif ((strfind(fileName, 'Prob3') == 1))
                obj.tempObj = Prob3Related(fileName);
            end
            
            % Initialize variables.
            obj.N_real_a = getRealN(obj.N_a);
            obj.N_real_b = getRealN(obj.N_b);
            obj.N_real_c = getRealN(obj.N_c);
            obj.N_real_d = getRealN(obj.N_d);
            
            mainCalculation(obj, mode);
        end
        
        % Get temperature at each node.
        % Matlab syntax: a function (other than the constructor) takes at
        % least one input, and the first input has to be the current class
        % (object) itself.
        function uxx = getUxx(obj)
            uxx = zeros(obj.totalNodes, 1);
            for i = 1 : obj.totalNodes
                uxx(i, 1) = obj.u_n(i * 2 - 1, 1);
            end
            return
        end
        
        function uyy = getUyy(obj)
            uyy = zeros(obj.totalNodes, 1);
            for i = 1 : obj.totalNodes
                uyy(i, 1) = obj.u_n(i * 2, 1);
            end
            return
        end
        
        function x = getX(obj)
            x = obj.nodeGlobal(:, 2);
            return
        end
        
        function y = getY(obj)
            y = obj.nodeGlobal(:, 3);
            return
        end
        
        % Get global nodes.
        function matrixN = getNodes(obj)
            matrixN = obj.nodeGlobal;
            return
        end
        
        % Get stress.
        function q = getStress(obj)
            q = obj.stress;
            return
        end
        
        % Get run time.
        function t = getTime(obj)
            t = obj.time;
            return
        end
        
        % Get elements.
        function ele = getElements(obj)
            ele = obj.elements;
            return
        end
        
        % Get relative residual array.
        function r = getRR(obj)
            r = obj.relativeResidual;
            return
        end
        
    end
    
    methods (Access = private)
        
        % Main calculation (private).
        function mainCalculation(obj, mode)
            
            % Start timing.
            %tic
            
            % Set the sizes of matrices.
            obj.K_global = zeros(obj.totalNodes * 2, obj.totalNodes * 2);
            obj.f_n = zeros(obj.totalNodes * 2, 1);
            obj.f_s = zeros(obj.totalNodes * 2, 1);
            obj.stress = zeros(obj.totalElements * 4, 6);
            obj.u_n = zeros(obj.totalNodes * 2, 1);
            obj.f_g = zeros(obj.totalNodes * 2, 1);
            
            % Initial guess for u_n.
            % Apply essential BC.
            for i = 1 : obj.tempObj.getEBCNum()
                currentEBC = obj.tempObj.getEssenBC(i); % The nodes.
                currentEssen = obj.tempObj.getEssen(i); % The nodal values.
                for j = 1 : size(currentEBC, 2)
                    tempX = obj.nodeGlobal(currentEBC(j), 2);
                    tempY = obj.nodeGlobal(currentEBC(j), 3);
                    tempEssen = currentEssen(tempX, tempY);
                    if (tempEssen(1, 1) ~= inf)
                        obj.u_n(currentEBC(j) * 2 - 1, 1) = tempEssen(1, 1);
                    end
                    if (tempEssen(2, 1) ~= inf)
                    obj.u_n(currentEBC(j) * 2, 1) = tempEssen(2, 1);
                    end
                end
            end
            
            % Iterations
            for n = 1 : obj.nMax
                
                % Zero out K_global, f_n, f_s, and f_g
                obj.K_global = zeros(obj.totalNodes * 2, obj.totalNodes * 2);
                obj.f_s = zeros(obj.totalNodes * 2, 1);
                obj.f_g = zeros(obj.totalNodes * 2, 1);
                
                % Start the main loop to calculate temperature.
                for i = 1 : obj.totalElements
                    
                    % Get the x, y nodal values for this element.
                    % Store the values in xy matrix.
                    % Find the gather matrix.
                    xy = zeros(8, 2);
                    L = zeros(8 * 2, obj.totalNodes * 2);
                    for j = 1 : 8
                        globalNodeNumber = obj.elements(i, j + 1);
                        xy(j, 1) = obj.nodeGlobal(globalNodeNumber, 2);
                        xy(j, 2) = obj.nodeGlobal(globalNodeNumber, 3);
                        L(j * 2 - 1, globalNodeNumber * 2 - 1) = 1;
                        L(j * 2, globalNodeNumber * 2) = 1;
                    end
                    
                    % Find K_global
                    % f_g is also calculated here implicitly.
                    obj.K_global = obj.K_global + getStiffnessMatrix(obj, xy, L, mode);
                    
                    % Apply Neumann BC.
                    if (n == 1)
                        obj.f_n = obj.f_n + getNeumannBC(obj, i, xy, L);
                    end
                    
                    % Apply body source.
                    obj.f_s = obj.f_s + getSource(obj, xy, L);
                end
                
                tempfns = obj.f_n + obj.f_s;
                                
                % Apply essential BC. The nodal values are always zeros for
                % those cases. This is a modified penalty method for those
                % specific cases.
                cMax = max(max(abs(obj.K_global)));
                cMax = cMax * 10 ^ 7;
                for k = 1 : obj.tempObj.getEBCNum()
                    currentEBC = obj.tempObj.getEssenBC(k); % The nodes.
                    currentEssen = obj.tempObj.getEssen(k); % The nodal values.
                    for j = 1 : size(currentEBC, 2)
                        tempX = obj.nodeGlobal(currentEBC(j), 2);
                        tempY = obj.nodeGlobal(currentEBC(j), 3);
                        tempEssen = currentEssen(tempX, tempY);
                        if (tempEssen(1, 1) ~= inf)
                            obj.K_global(currentEBC(j) * 2 - 1, currentEBC(j) * 2 - 1) = cMax;
                            tempfns(currentEBC(j) * 2 - 1, 1) = obj.f_g(currentEBC(j) * 2 - 1, 1);
                        end
                        if (tempEssen(2, 1) ~= inf)
                            obj.K_global(currentEBC(j) * 2, currentEBC(j) * 2) = cMax;
                            tempfns(currentEBC(j) * 2, 1) = obj.f_g(currentEBC(j) * 2, 1);
                        end
                    end
                end
                
                % Combine the RHS.
                RHS = tempfns - obj.f_g;
                tempRe = norm(RHS) / norm(tempfns);
                
                if (tempRe < obj.err)
                    break
                end
                tempU = obj.K_global \ RHS;
                obj.u_n = obj.u_n + tempU;
                obj.relativeResidual = [obj.relativeResidual, tempRe];
            end
            
            if (n == obj.nMax)
                convergence = 0
            end
            
            % Start the main loop to calculate stress at quadrature points.
            for i = 1 : obj.totalElements
                
                % Get the x, y nodal values for this element.
                % Store the values in xy matrix.
                % Find the temperature at local nodes.
                xy = zeros(8, 2);
                uxuy = zeros(8, 2);
                for j = 1 : 8
                    globalNodeNumber = obj.elements(i, j + 1);
                    xy(j, 1) = obj.nodeGlobal(globalNodeNumber, 2);
                    xy(j, 2) = obj.nodeGlobal(globalNodeNumber, 3);
                    uxuy(j, 1) = obj.u_n(globalNodeNumber * 2 - 1, 1);
                    uxuy(j, 2) = obj.u_n(globalNodeNumber * 2, 1);
                end
                
                % Iterate through 4 quadrature points in each elements.
                stressLocal = localStressNonLinear(xy, obj.N_a, obj.N_b, obj.N_c, obj.N_d, obj.G_a, obj.G_b, obj.G_c, obj.G_d, uxuy, obj.v0, obj.E0, obj.q0);
                for j = 1 : 4
                    obj.stress(i * 4 - 3 + j - 1, :) = stressLocal(j, :);
                end
            end

%             % End timing
%             %obj.time = toc;
        end
        
        % Get stifness matrix.
        % Also get gForce matrix here.
        % Private.
        function K_el = getStiffnessMatrix(obj, xy, L, mode)
            
            % Find the jacobian matrices and their determinants.
            % 2 x 2
            J_a = obj.G_a * xy;
            detJ_a = det(J_a);
            J_b = obj.G_b * xy;
            detJ_b = det(J_b);
            J_c = obj.G_c * xy;
            detJ_c = det(J_c);
            J_d = obj.G_d * xy;
            detJ_d = det(J_d);
            
            % Find B matrices. 2 x 8
            B_a = inv(J_a) * obj.G_a;
            B_b = inv(J_b) * obj.G_b;
            B_c = inv(J_c) * obj.G_c;
            B_d = inv(J_d) * obj.G_d;
            
            % Real B matrices. 3 x 16
            B_real_a = getRealB(B_a);
            B_real_b = getRealB(B_b);
            B_real_c = getRealB(B_c);
            B_real_d = getRealB(B_d);
            
            % Epsilon matrix: eps_xx, eps_yy, gamma_xy
            eps_a = B_real_a * L * obj.u_n;
            eps_b = B_real_b * L * obj.u_n;
            eps_c = B_real_c * L * obj.u_n;
            eps_d = B_real_d * L * obj.u_n;
            
            % eij matrix. 2 x 2.
            eij_a = zeros(2, 2);
            eij_a(1, 1) = eps_a(1, 1) - 1 / 3 * (eps_a(1, 1) + eps_a(2, 1));
            eij_a(2, 2) = eps_a(2, 1) - 1 / 3 * (eps_a(1, 1) + eps_a(2, 1));
            eij_a(1, 2) = eps_a(3, 1);
            eij_a(2, 1) = eij_a(1, 2);
            
            eij_b = zeros(2, 2);
            eij_b(1, 1) = eps_b(1, 1) - 1 / 3 * (eps_b(1, 1) + eps_b(2, 1));
            eij_b(2, 2) = eps_b(2, 1) - 1 / 3 * (eps_b(1, 1) + eps_b(2, 1));
            eij_b(1, 2) = eps_b(3, 1);
            eij_b(2, 1) = eij_b(1, 2);
            
            eij_c = zeros(2, 2);
            eij_c(1, 1) = eps_c(1, 1) - 1 / 3 * (eps_c(1, 1) + eps_c(2, 1));
            eij_c(2, 2) = eps_c(2, 1) - 1 / 3 * (eps_c(1, 1) + eps_c(2, 1));
            eij_c(1, 2) = eps_c(3, 1);
            eij_c(2, 1) = eij_c(1, 2);
            
            eij_d = zeros(2, 2);
            eij_d(1, 1) = eps_d(1, 1) - 1 / 3 * (eps_d(1, 1) + eps_d(2, 1));
            eij_d(2, 2) = eps_d(2, 1) - 1 / 3 * (eps_d(1, 1) + eps_d(2, 1));
            eij_d(1, 2) = eps_d(3, 1);
            eij_d(2, 1) = eij_d(1, 2);
            
            % Scalar.
            epse_a = (2 / 3 * ((eij_a(1, 1)) ^ 2 + (eij_a(2, 2)) ^ 2 + (eij_a(1, 2)) ^ 2 + (eij_a(2, 1)) ^ 2)) ^ 0.5;
            epse_b = (2 / 3 * ((eij_b(1, 1)) ^ 2 + (eij_b(2, 2)) ^ 2 + (eij_b(1, 2)) ^ 2 + (eij_b(2, 1)) ^ 2)) ^ 0.5;
            epse_c = (2 / 3 * ((eij_c(1, 1)) ^ 2 + (eij_c(2, 2)) ^ 2 + (eij_c(1, 2)) ^ 2 + (eij_c(2, 1)) ^ 2)) ^ 0.5;
            epse_d = (2 / 3 * ((eij_d(1, 1)) ^ 2 + (eij_d(2, 2)) ^ 2 + (eij_d(1, 2)) ^ 2 + (eij_d(2, 1)) ^ 2)) ^ 0.5;
            
            if (mode == 0) % tangent
                kT_a = obj.E0 * (1 + epse_a ^ obj.q0) / 3 / (1 - 2 * obj.v0);
                GT_a = obj.E0 * (1 + (obj.q0 + 1) * epse_a ^ obj.q0) / 2 / (1 + obj.v0);
                vD_a = (3 * kT_a - 2 * GT_a) / (3 * kT_a + GT_a) / 2;
                ED_a = 2 * GT_a * (1 + vD_a);
                
                kT_b = obj.E0 * (1 + epse_b ^ obj.q0) / 3 / (1 - 2 * obj.v0);
                GT_b = obj.E0 * (1 + (obj.q0 + 1) * epse_b ^ obj.q0) / 2 / (1 + obj.v0);
                vD_b = (3 * kT_b - 2 * GT_b) / (3 * kT_b + GT_b) / 2;
                ED_b = 2 * GT_b * (1 + vD_b);
                
                kT_c = obj.E0 * (1 + epse_c ^ obj.q0) / 3 / (1 - 2 * obj.v0);
                GT_c = obj.E0 * (1 + (obj.q0 + 1) * epse_c ^ obj.q0) / 2 / (1 + obj.v0);
                vD_c = (3 * kT_c - 2 * GT_c) / (3 * kT_c + GT_c) / 2;
                ED_c = 2 * GT_c * (1 + vD_c);
                
                kT_d = obj.E0 * (1 + epse_d ^ obj.q0) / 3 / (1 - 2 * obj.v0);
                GT_d = obj.E0 * (1 + (obj.q0 + 1) * epse_d ^ obj.q0) / 2 / (1 + obj.v0);
                vD_d = (3 * kT_d - 2 * GT_d) / (3 * kT_d + GT_d) / 2;
                ED_d = 2 * GT_d * (1 + vD_d);
                
                % D matrix
                D_a = getDMatrix(vD_a, ED_a);
                D_b = getDMatrix(vD_b, ED_b);
                D_c = getDMatrix(vD_c, ED_c);
                D_d = getDMatrix(vD_d, ED_d);
                
                vDs_a = obj.v0;
                EDs_a = obj.E0 * (1 + epse_a ^ obj.q0);
                
                vDs_b = obj.v0;
                EDs_b = obj.E0 * (1 + epse_b ^ obj.q0);
                
                vDs_c = obj.v0;
                EDs_c = obj.E0 * (1 + epse_c ^ obj.q0);
                
                vDs_d = obj.v0;
                EDs_d = obj.E0 * (1 + epse_d ^ obj.q0);
                
                % Dr matrix
                Dr_a = getDMatrix(vDs_a, EDs_a);
                Dr_b = getDMatrix(vDs_b, EDs_b);
                Dr_c = getDMatrix(vDs_c, EDs_c);
                Dr_d = getDMatrix(vDs_d, EDs_d);                
            else % secant
                vD_a = obj.v0;
                ED_a = obj.E0 * (1 + epse_a ^ obj.q0);
                
                vD_b = obj.v0;
                ED_b = obj.E0 * (1 + epse_b ^ obj.q0);
                
                vD_c = obj.v0;
                ED_c = obj.E0 * (1 + epse_c ^ obj.q0);
                
                vD_d = obj.v0;
                ED_d = obj.E0 * (1 + epse_d ^ obj.q0);
                
                % D matrix
                D_a = getDMatrix(vD_a, ED_a);
                D_b = getDMatrix(vD_b, ED_b);
                D_c = getDMatrix(vD_c, ED_c);
                D_d = getDMatrix(vD_d, ED_d);
                
                % Dr matrix
                Dr_a = D_a;
                Dr_b = D_b;
                Dr_c = D_c;
                Dr_d = D_d;
            end
            
            K_el = B_real_a.' * D_a * B_real_a * detJ_a + B_real_b.' * D_b * B_real_b * detJ_b + B_real_c.' * D_c * B_real_c * detJ_c + B_real_d.' * D_d * B_real_d * detJ_d;
            temp_fg = B_real_a.' * Dr_a * eps_a * detJ_a + B_real_b.' * Dr_b * eps_b * detJ_b + B_real_c.' * Dr_c * eps_c * detJ_c + B_real_d.' * Dr_d * eps_d * detJ_d;
            
            % Transfer local nodes into global nodes.
            K_el = L.' * K_el * L;
            temp_fg = L.' * temp_fg;
            obj.f_g = obj.f_g + temp_fg;
            return
        end
        
        % Get Neumann BC. (private)
        function force_n = getNeumannBC(obj, elementNum, xy, L)
            force_n = zeros(obj.totalNodes * 2, 1);
            currentElement = obj.elements(elementNum, :);
            
            % Go through each Neumann BC for each element.
            for i = 1 : obj.tempObj.getnBCNum()
                
                % Check if a node is on the boundary.
                % Use a 1x2 matrix to store the node number.
                % Each element can only have 0, 1, or 2 corner nodes on a
                % certain boundary.
                % loop after k reaches 3.
                currentBoundary = obj.tempObj.getNBC(i);
                currentFlux = obj.tempObj.getFlux(i);
                record = zeros(1, 2);
                k = 1;
                for j = 1 : 4
                    currentElement(1, j + 1);%
                    if ismember(currentElement(j + 1), currentBoundary) == 1
                        record(1, k) = j;
                        k = k + 1;
                    end
                    if k == 3
                        break
                    end
                end
                
                recordState = record(1) * record(2);
                if recordState ~= 0
                    force_n = force_n + force2PointGaussQ8(xy, L, recordState, currentFlux);
                end
            end
            %force_n = -1 * force_n; % !!!!! here not above. Don't need it.
            return
        end
        
        % Get body source.
        % private.
        function force_b = getSource(obj, xy, L)
            force_b = zeros(obj.totalNodes * 2, 1);
            J_a = obj.G_a * xy;
            detJ_a = det(J_a);
            J_b = obj.G_b * xy;
            detJ_b = det(J_b);
            J_c = obj.G_c * xy;
            detJ_c = det(J_c);
            J_d = obj.G_d * xy;
            detJ_d = det(J_d);
            
            % Go through each source for each element.
            for i = 1 : obj.tempObj.getSourceNum()
                currentSource = obj.tempObj.getSource(i);
                force_b = force_b + source2PointGaussQ8(obj.N_a, obj.N_b, obj.N_c, obj.N_d, detJ_a, detJ_b, detJ_c, detJ_d, L, xy, currentSource, obj.N_real_a, obj.N_real_b, obj.N_real_c, obj.N_real_d);
            end
            return
        end
        
    end
end
