classdef UniformTPM
    properties
        U
        L
    end
    methods
        function obj = UniformTPM(L,U)
            if nargin == 2
                obj.U = U;
                obj.L = L;
            end
        end
        function e = X(obj)
            e = (obj.U+obj.L)/2;
        end
        function e = XX(obj)
            % E[x^2]
            e = (obj.U^2+obj.U*obj.L+obj.L^2)/3;
        end
        function e = C(obj)
            % E[cosx]
            e = (sin(obj.U)-sin(obj.L))/(obj.U-obj.L);
        end
        function e = S(obj)
            % E[sinx]
            e = -(cos(obj.U)-cos(obj.L))/(obj.U-obj.L);
        end
        function e = XC(obj)
            % E[xcosx]
            e = -(-obj.U*sin(obj.U)+obj.L*sin(obj.L)-cos(obj.U)+cos(obj.L))/(obj.U-obj.L);
        end
        function e = XS(obj)
            % E[xsinx]
            e = -(obj.U*cos(obj.U)-obj.L*cos(obj.L)-sin(obj.U)+sin(obj.L))/(obj.U-obj.L);
        end
        function e = XXC(obj)
            % E[x^2 cosx]
            e = -((-obj.U^2+2)*sin(obj.U)-2*obj.U*cos(obj.U)+(obj.L^2-2)*sin(obj.L)+2*obj.L*cos(obj.L))/(obj.U-obj.L);
        end
        function e = XXS(obj)
            % E[x^2 sinx]
            e =((-obj.U^2+2)*cos(obj.U)+2*obj.U*sin(obj.U)+(obj.L^2-2)*cos(obj.L)-2*obj.L*sin(obj.L))/(obj.U-obj.L);
        end
        function e = CS(obj)
            e = -0.5*(cos(2*obj.U)-cos(2*obj.L))/(2*obj.U-2*obj.L);
        end
        function e = CC(obj)
            e = 0.5*((sin(2*obj.U)-sin(2*obj.L))/(2*obj.U-2*obj.L)+1);
        end
        function e = SS(obj)
            e = 0.5*(1-(sin(2*obj.U)-sin(2*obj.L))/(2*obj.U-2*obj.L));
        end
    end
end
