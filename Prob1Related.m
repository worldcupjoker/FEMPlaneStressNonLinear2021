% Copyright, 2021, (Ted) Zhongkun Zhan, all rights reserved.
% File name: Prob1Q8Related

classdef Prob1Related < handle
    properties (Access = private)
        nBCNum = 1;
        essenBCNum = 2;
        sourceNum = 1;
        
        nBC1;
        nBC2;
        essenBC1; % the nodes indices.
        essenBC2;
        essen1 = [0; Inf]; % 0 is the real BC; Inf means not BC.
        essen2 = [Inf; 0];
        flux1 = [20e6; 0];
        bodySource = [0; 0];
    end
    
    methods
        
        % Constructor
        function obj = Prob1Related(fileName)
            obj.nBC1 = readinp('*Nset, nset=Right', fileName);
            obj.essenBC1 = readinp('*Nset, nset=Left', fileName);
            obj.essenBC2 = readinp('*Nset, nset=Bottom', fileName);
        end
        
        function neumannBC = getNBC(obj, i)
            if i == 1
                neumannBC = obj.nBC1;
                return
            end
            if i == 2
                neumannBC = obj.nBC2;
                return
            end
        end
        
        function flux = getFlux(obj, i)
            if i == 1
                flux = @(x, y)(obj.flux1);
                return
            end
            if i == 2
                flux = @(x, y)(obj.flux1);
                return
            end
        end
        
        function s = getSource(obj, i)
            if i == 1
                s = @(x, y)(obj.bodySource);
                return
            end
        end
        
        function tBC = getEssenBC(obj, i)
            if i == 1
                tBC = obj.essenBC1;
                return
            end
            if i == 2
                tBC = obj.essenBC2;
                return
            end
        end
        
        function temp = getEssen(obj, i)
            if i == 1
                temp = @(x, y)(obj.essen1);
                return
            end
            if i == 2
                temp = @(x, y)(obj.essen2);
                return
            end
        end
        
        function num = getnBCNum(obj)
            num = obj.nBCNum;
            return
        end
        
        function num = getSourceNum(obj)
            num = obj.sourceNum;
            return
        end
        
        function num = getEBCNum(obj)
            num = obj.essenBCNum;
            return
        end
        
    end
end