classdef ExpotensialTPM
    properties
        lambda        
    end
    methods
        function obj = ExpotensialTPM(lambda)
            if nargin == 1
                obj.lambda = lambda;
            end
        end
        function e = X(obj)
            e = 1/obj.lambda;
        end
        function e = XX(obj)
            % E[x^2]
            e = 2/obj.lambda^2;
        end
        function e = C(obj)
            % E[cosx]
            e = obj.lambda^2/(obj.lambda^2+1);
        end
        function e = S(obj)
            % E[sinx]
            e = obj.lambda/(obj.lambda^2+1);
        end
        function e = XC(obj)
            % E[xcosx]
            e = (obj.lambda^3-obj.lambda)/(obj.lambda^2+1)^2;
        end
        function e = XS(obj)
            % E[xsinx]
            e = (2*obj.lambda^2)/(obj.lambda^2+1)^2;
        end
        function e = XXC(obj)
            % E[x^2 cosx]
            e = (2*obj.lambda^2)*(obj.lambda^2-3)/(obj.lambda^2+1)^3;
        end
        function e = XXS(obj)
            % E[x^2 sinx]
            e =(2*obj.lambda)*(3*obj.lambda^2-1)/(obj.lambda^2+1)^3;
        end
    end
end
