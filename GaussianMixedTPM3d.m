classdef GaussianMixedTPM3d < GaussianTPM
    properties
        Mu,Sigma2
        T,Lambda
        Z31,Z32,Z33
    end
    methods
        function obj = GaussianMixedTPM3d(Mu,Sigma2)
            if nargin == 2
                if size(Mu,1) ~= 3
                    error('Not 3-dimension, check dimension of mu or sigma')
                end
                obj.Mu = Mu;
                obj.Sigma2 = Sigma2;
                [obj.T,obj.Lambda] = eig(Sigma2);
                T31 = obj.T(3,1);
                T32 = obj.T(3,2);
                T33 = obj.T(3,3);
                mu_y = obj.T'*obj.Mu;
                mu_y1 = mu_y(1);
                mu_y2 = mu_y(2);
                mu_y3 = mu_y(3);
                lambda1 = obj.Lambda(1,1);
                lambda2 = obj.Lambda(2,2);
                lambda3 = obj.Lambda(3,3);
                % 確率分布のクラス定義
                obj.Z31 = GaussianTPM(T31*mu_y1,(T31^2)*lambda1);
                obj.Z32 = GaussianTPM(T32*mu_y2,(T32^2)*lambda2);
                obj.Z33 = GaussianTPM(T33*mu_y3,(T33^2)*lambda3);
            end
        end
        function e = X(obj)
            e = obj.Mu(1);
        end
        function e = Y(obj)
            e = obj.Mu(2);
        end
        function e = Th(obj)
            e = obj.Mu(3);
        end

        function e = XC(obj)
            % calc E[(x*cos(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(1)...
                    *(cos(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)*(z31.XCosX*z32.CosX*z33.CosX-z31.XCosX*z32.SinX*z33.SinX-z31.XSinX*z32.CosX*z33.SinX-z31.XSinX*z32.SinX*z33.CosX)...
                    +(T12/T32)*(z31.CosX*z32.XCosX*z33.CosX-z31.CosX*z32.XSinX*z33.SinX-z31.SinX*z32.XCosX*z33.SinX-z31.SinX*z32.XSinX*z33.CosX)...
                    +(T13/T33)*(z31.CosX*z32.CosX*z33.XCosX-z31.CosX*z32.SinX*z33.XSinX-z31.SinX*z32.CosX*z33.XSinX-z31.SinX*z32.SinX*z33.XCosX);
            end
        end
        function e = XS(obj)
            % calc E[(x*sin(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(1)...
                    *(sin(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)*(-z31.XSinX*z32.SinX*z33.SinX+z31.XSinX*z32.CosX*z33.CosX+z31.XCosX*z32.SinX*z33.CosX+z31.XCosX*z32.CosX*z33.SinX)...
                    +(T12/T32)*(-z31.SinX*z32.XSinX*z33.SinX+z31.SinX*z32.XCosX*z33.CosX+z31.CosX*z32.XSinX*z33.CosX+z31.CosX*z32.XCosX*z33.SinX)...
                    +(T13/T33)*(-z31.SinX*z32.SinX*z33.XSinX+z31.SinX*z32.CosX*z33.XCosX+z31.CosX*z32.SinX*z33.XCosX+z31.CosX*z32.CosX*z33.XSinX);
            end
        end
        function e = YC(obj)
            % calc E[(y*cos(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(2)...
                    *(cos(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)*(z31.XCosX*z32.CosX*z33.CosX-z31.XCosX*z32.SinX*z33.SinX-z31.XSinX*z32.CosX*z33.SinX-z31.XSinX*z32.SinX*z33.CosX)...
                    +(T22/T32)*(z31.CosX*z32.XCosX*z33.CosX-z31.CosX*z32.XSinX*z33.SinX-z31.SinX*z32.XCosX*z33.SinX-z31.SinX*z32.XSinX*z33.CosX)...
                    +(T23/T33)*(z31.CosX*z32.CosX*z33.XCosX-z31.CosX*z32.SinX*z33.XSinX-z31.SinX*z32.CosX*z33.XSinX-z31.SinX*z32.SinX*z33.XCosX);
            end
        end
        function e = YS(obj)
            % calc E[(y*sin(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(2)...
                    *(sin(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)*(-z31.XSinX*z32.SinX*z33.SinX+z31.XSinX*z32.CosX*z33.CosX+z31.XCosX*z32.SinX*z33.CosX+z31.XCosX*z32.CosX*z33.SinX)...
                    +(T22/T32)*(-z31.SinX*z32.XSinX*z33.SinX+z31.SinX*z32.XCosX*z33.CosX+z31.CosX*z32.XSinX*z33.CosX+z31.CosX*z32.XCosX*z33.SinX)...
                    +(T23/T33)*(-z31.SinX*z32.SinX*z33.XSinX+z31.SinX*z32.CosX*z33.XCosX+z31.CosX*z32.SinX*z33.XCosX+z31.CosX*z32.CosX*z33.XSinX);
            end
        end
        function e = X2C(obj)
            % calc E[(x*x*cos(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = (obj.Sigma2(1,1)+obj.Mu(1)^2)...
                    *(cos(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)^2*(z31.X2CosX*z32.CosX*z33.CosX-z31.X2CosX*z32.SinX*z33.SinX-z31.X2SinX*z32.CosX*z33.SinX-z31.X2SinX*z32.SinX*z33.CosX)...
                    +(T12/T32)^2*(z31.CosX*z32.X2CosX*z33.CosX-z31.CosX*z32.X2SinX*z33.SinX-z31.SinX*z32.X2CosX*z33.SinX-z31.SinX*z32.X2SinX*z33.CosX)...
                    +(T13/T33)^2*(z31.CosX*z32.CosX*z33.X2CosX-z31.CosX*z32.SinX*z33.X2SinX-z31.SinX*z32.CosX*z33.X2SinX-z31.SinX*z32.SinX*z33.X2CosX)...
                    +2*(T11/T31)*(T12/T32)*(z31.XCosX*z32.XCosX*z33.CosX-z31.XCosX*z32.XSinX*z33.SinX-z31.XSinX*z32.XCosX*z33.SinX-z31.XSinX*z32.XSinX*z33.CosX)...
                    +2*(T12/T32)*(T13/T33)*(z31.CosX*z32.XCosX*z33.XCosX-z31.CosX*z32.XSinX*z33.XSinX-z31.SinX*z32.XCosX*z33.XSinX-z31.SinX*z32.XSinX*z33.XCosX)...
                    +2*(T13/T33)*(T11/T31)*(z31.XCosX*z32.CosX*z33.XCosX-z31.XCosX*z32.SinX*z33.XSinX-z31.XSinX*z32.CosX*z33.XSinX-z31.XSinX*z32.SinX*z33.XCosX);
            end
        end
        function e = X2S(obj)
            % calc E[(x*x*sin(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = (obj.Sigma2(1,1)+obj.Mu(1)^2)...
                    *(sin(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)^2*(-z31.X2SinX*z32.SinX*z33.SinX+z31.X2SinX*z32.CosX*z33.CosX+z31.X2CosX*z32.SinX*z33.CosX+z31.X2CosX*z32.CosX*z33.SinX)...
                    +(T12/T32)^2*(-z31.SinX*z32.X2SinX*z33.SinX+z31.SinX*z32.X2CosX*z33.CosX+z31.CosX*z32.X2SinX*z33.CosX+z31.CosX*z32.X2CosX*z33.SinX)...
                    +(T13/T33)^2*(-z31.SinX*z32.SinX*z33.X2SinX+z31.SinX*z32.CosX*z33.X2CosX+z31.CosX*z32.SinX*z33.X2CosX+z31.CosX*z32.CosX*z33.X2SinX)...
                    +2*(T11/T31)*(T12/T32)*(-z31.XSinX*z32.XSinX*z33.SinX+z31.XSinX*z32.XCosX*z33.CosX+z31.XCosX*z32.XSinX*z33.CosX+z31.XCosX*z32.XCosX*z33.SinX)...
                    +2*(T12/T32)*(T13/T33)*(-z31.SinX*z32.XSinX*z33.XSinX+z31.SinX*z32.XCosX*z33.XCosX+z31.CosX*z32.XSinX*z33.XCosX+z31.CosX*z32.XCosX*z33.XSinX)...
                    +2*(T13/T33)*(T11/T31)*(-z31.XSinX*z32.SinX*z33.XSinX+z31.XSinX*z32.CosX*z33.XCosX+z31.XCosX*z32.SinX*z33.XCosX+z31.XCosX*z32.CosX*z33.XSinX);
            end
        end
        function e = Y2C(obj)
            % calc E[(y*y*cos(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = (obj.Sigma2(2,2)+obj.Mu(2)^2)...
                    *(cos(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)^2*(z31.X2CosX*z32.CosX*z33.CosX-z31.X2CosX*z32.SinX*z33.SinX-z31.X2SinX*z32.CosX*z33.SinX-z31.X2SinX*z32.SinX*z33.CosX)...
                    +(T22/T32)^2*(z31.CosX*z32.X2CosX*z33.CosX-z31.CosX*z32.X2SinX*z33.SinX-z31.SinX*z32.X2CosX*z33.SinX-z31.SinX*z32.X2SinX*z33.CosX)...
                    +(T23/T33)^2*(z31.CosX*z32.CosX*z33.X2CosX-z31.CosX*z32.SinX*z33.X2SinX-z31.SinX*z32.CosX*z33.X2SinX-z31.SinX*z32.SinX*z33.X2CosX)...
                    +2*(T21/T31)*(T22/T32)*(z31.XCosX*z32.XCosX*z33.CosX-z31.XCosX*z32.XSinX*z33.SinX-z31.XSinX*z32.XCosX*z33.SinX-z31.XSinX*z32.XSinX*z33.CosX)...
                    +2*(T22/T32)*(T23/T33)*(z31.CosX*z32.XCosX*z33.XCosX-z31.CosX*z32.XSinX*z33.XSinX-z31.SinX*z32.XCosX*z33.XSinX-z31.SinX*z32.XSinX*z33.XCosX)...
                    +2*(T23/T33)*(T21/T31)*(z31.XCosX*z32.CosX*z33.XCosX-z31.XCosX*z32.SinX*z33.XSinX-z31.XSinX*z32.CosX*z33.XSinX-z31.XSinX*z32.SinX*z33.XCosX);
            end
        end
        function e = Y2S(obj)
            % calc E[(y*y*sin(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = (obj.Sigma2(2,2)+obj.Mu(2)^2)...
                    *(sin(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)^2*(-z31.X2SinX*z32.SinX*z33.SinX+z31.X2SinX*z32.CosX*z33.CosX+z31.X2CosX*z32.SinX*z33.CosX+z31.X2CosX*z32.CosX*z33.SinX)...
                    +(T22/T32)^2*(-z31.SinX*z32.X2SinX*z33.SinX+z31.SinX*z32.X2CosX*z33.CosX+z31.CosX*z32.X2SinX*z33.CosX+z31.CosX*z32.X2CosX*z33.SinX)...
                    +(T23/T33)^2*(-z31.SinX*z32.SinX*z33.X2SinX+z31.SinX*z32.CosX*z33.X2CosX+z31.CosX*z32.SinX*z33.X2CosX+z31.CosX*z32.CosX*z33.X2SinX)...
                    +2*(T21/T31)*(T22/T32)*(-z31.XSinX*z32.XSinX*z33.SinX+z31.XSinX*z32.XCosX*z33.CosX+z31.XCosX*z32.XSinX*z33.CosX+z31.XCosX*z32.XCosX*z33.SinX)...
                    +2*(T22/T32)*(T23/T33)*(-z31.SinX*z32.XSinX*z33.XSinX+z31.SinX*z32.XCosX*z33.XCosX+z31.CosX*z32.XSinX*z33.XCosX+z31.CosX*z32.XCosX*z33.XSinX)...
                    +2*(T23/T33)*(T21/T31)*(-z31.XSinX*z32.SinX*z33.XSinX+z31.XSinX*z32.CosX*z33.XCosX+z31.XCosX*z32.SinX*z33.XCosX+z31.XCosX*z32.CosX*z33.XSinX);
            end
        end
        function e = XYC(obj)
            % calc E[(x*y*cos(theta))]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = (obj.Mu(1))...
                    *(obj.Mu(2))...
                    *(cos(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(T21/T31)*(z31.X2CosX*z32.CosX*z33.CosX-z31.X2CosX*z32.SinX*z33.SinX-z31.X2SinX*z32.CosX*z33.SinX-z31.X2SinX*z32.SinX*z33.CosX)...
                    +(T12/T32)*(T22/T32)*(z31.CosX*z32.X2CosX*z33.CosX-z31.CosX*z32.X2SinX*z33.SinX-z31.SinX*z32.X2CosX*z33.SinX-z31.SinX*z32.X2SinX*z33.CosX)...
                    +(T13/T33)*(T23/T33)*(z31.CosX*z32.CosX*z33.X2CosX-z31.CosX*z32.SinX*z33.X2SinX-z31.SinX*z32.CosX*z33.X2SinX-z31.SinX*z32.SinX*z33.X2CosX)...
                    +((T11/T31)*(T22/T32)+(T12/T32)*(T21/T31))*(z31.XCosX*z32.XCosX*z33.CosX-z31.XCosX*z32.XSinX*z33.SinX-z31.XSinX*z32.XCosX*z33.SinX-z31.XSinX*z32.XSinX*z33.CosX)...
                    +((T12/T32)*(T23/T33)+(T13/T33)*(T22/T32))*(z31.CosX*z32.XCosX*z33.XCosX-z31.CosX*z32.XSinX*z33.XSinX-z31.SinX*z32.XCosX*z33.XSinX-z31.SinX*z32.XSinX*z33.XCosX)...
                    +((T13/T33)*(T21/T31)+(T11/T31)*(T23/T33))*(z31.XCosX*z32.CosX*z33.XCosX-z31.XCosX*z32.SinX*z33.XSinX-z31.XSinX*z32.CosX*z33.XSinX-z31.XSinX*z32.SinX*z33.XCosX);
            end
        end
        function e = XYS(obj)
            % calc E[(x*y*sin(theta))]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = (obj.Mu(1))...
                    *(obj.Mu(2))...
                    *(sin(obj.Mu(3))*exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(T21/T31)*(-z31.X2SinX*z32.SinX*z33.SinX+z31.X2SinX*z32.CosX*z33.CosX+z31.X2CosX*z32.SinX*z33.CosX+z31.X2CosX*z32.CosX*z33.SinX)...
                    +(T12/T32)*(T22/T32)*(-z31.SinX*z32.X2SinX*z33.SinX+z31.SinX*z32.X2CosX*z33.CosX+z31.CosX*z32.X2SinX*z33.CosX+z31.CosX*z32.X2CosX*z33.SinX)...
                    +(T13/T33)*(T23/T33)*(-z31.SinX*z32.SinX*z33.X2SinX+z31.SinX*z32.CosX*z33.X2CosX+z31.CosX*z32.SinX*z33.X2CosX+z31.CosX*z32.CosX*z33.X2SinX)...
                    +((T11/T31)*(T22/T32)+(T12/T32)*(T21/T31))*(-z31.XSinX*z32.XSinX*z33.SinX+z31.XSinX*z32.XCosX*z33.CosX+z31.XCosX*z32.XSinX*z33.CosX+z31.XCosX*z32.XCosX*z33.SinX)...
                    +((T12/T32)*(T23/T33)+(T13/T33)*(T22/T32))*(-z31.SinX*z32.XSinX*z33.XSinX+z31.SinX*z32.XCosX*z33.XCosX+z31.CosX*z32.XSinX*z33.XCosX+z31.CosX*z32.XCosX*z33.XSinX)...
                    +((T13/T33)*(T21/T31)+(T11/T31)*(T23/T33))*(-z31.XSinX*z32.SinX*z33.XSinX+z31.XSinX*z32.CosX*z33.XCosX+z31.XCosX*z32.SinX*z33.XCosX+z31.XCosX*z32.CosX*z33.XSinX);
            end
        end
        function e=XThC(obj)
            % calc E[(x*theta*cos(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(1)...
                    *((obj.Mu(3)*cos(obj.Mu(3))-obj.Sigma2(3,3)*sin(obj.Mu(3))) * exp(-0.5*obj.Sigma2(3,3)));
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(z31.X2CosX*z32.CosX*z33.CosX-z31.X2CosX*z32.SinX*z33.SinX-z31.X2SinX*z32.CosX*z33.SinX-z31.X2SinX*z32.SinX*z33.CosX)...
                    +(T12/T32)*(z31.CosX*z32.X2CosX*z33.CosX-z31.CosX*z32.X2SinX*z33.SinX-z31.SinX*z32.X2CosX*z33.SinX-z31.SinX*z32.X2SinX*z33.CosX)...
                    +(T13/T33)*(z31.CosX*z32.CosX*z33.X2CosX-z31.CosX*z32.SinX*z33.X2SinX-z31.SinX*z32.CosX*z33.X2SinX-z31.SinX*z32.SinX*z33.X2CosX)...
                    +((T11/T31)+(T12/T32))*(z31.XCosX*z32.XCosX*z33.CosX-z31.XCosX*z32.XSinX*z33.SinX-z31.XSinX*z32.XCosX*z33.SinX-z31.XSinX*z32.XSinX*z33.CosX)...
                    +((T12/T32)+(T13/T33))*(z31.CosX*z32.XCosX*z33.XCosX-z31.CosX*z32.XSinX*z33.XSinX-z31.SinX*z32.XCosX*z33.XSinX-z31.SinX*z32.XSinX*z33.XCosX)...
                    +((T11/T31)+(T13/T33))*(z31.XCosX*z32.CosX*z33.XCosX-z31.XCosX*z32.SinX*z33.XSinX-z31.XSinX*z32.CosX*z33.XSinX-z31.XSinX*z32.SinX*z33.XCosX);
            end
        end
        function e=XThS(obj)
            % calc E[(x*theta*sin(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(1)...
                    *((obj.Mu(3)*sin(obj.Mu(3))+obj.Sigma2(3,3)*cos(obj.Mu(3))) * exp(-0.5*obj.Sigma2(3,3)));
            else

                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(-z31.X2SinX*z32.SinX*z33.SinX+z31.X2SinX*z32.CosX*z33.CosX+z31.X2CosX*z32.SinX*z33.CosX+z31.X2CosX*z32.CosX*z33.SinX)...
                    +(T12/T32)*(-z31.SinX*z32.X2SinX*z33.SinX+z31.SinX*z32.X2CosX*z33.CosX+z31.CosX*z32.X2SinX*z33.CosX+z31.CosX*z32.X2CosX*z33.SinX)...
                    +(T13/T33)*(-z31.SinX*z32.SinX*z33.X2SinX+z31.SinX*z32.CosX*z33.X2CosX+z31.CosX*z32.SinX*z33.X2CosX+z31.CosX*z32.CosX*z33.X2SinX)...
                    +((T11/T31)+(T12/T32))*(-z31.XSinX*z32.XSinX*z33.SinX+z31.XSinX*z32.XCosX*z33.CosX+z31.XCosX*z32.XSinX*z33.CosX+z31.XCosX*z32.XCosX*z33.SinX)...
                    +((T12/T32)+(T13/T33))*(-z31.SinX*z32.XSinX*z33.XSinX+z31.SinX*z32.XCosX*z33.XCosX+z31.CosX*z32.XSinX*z33.XCosX+z31.CosX*z32.XCosX*z33.XSinX)...
                    +((T11/T31)+(T13/T33))*(-z31.XSinX*z32.SinX*z33.XSinX+z31.XSinX*z32.CosX*z33.XCosX+z31.XCosX*z32.SinX*z33.XCosX+z31.XCosX*z32.CosX*z33.XSinX);
            end
        end
        function e=YThC(obj)
            % calc E[(y*theta*cos(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(2)...
                    *((obj.Mu(3)*cos(obj.Mu(3))-obj.Sigma2(3,3)*sin(obj.Mu(3))) * exp(-0.5*obj.Sigma2(3,3)));
            else

                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T21/T31)*(z31.X2CosX*z32.CosX*z33.CosX-z31.X2CosX*z32.SinX*z33.SinX-z31.X2SinX*z32.CosX*z33.SinX-z31.X2SinX*z32.SinX*z33.CosX)...
                    +(T22/T32)*(z31.CosX*z32.X2CosX*z33.CosX-z31.CosX*z32.X2SinX*z33.SinX-z31.SinX*z32.X2CosX*z33.SinX-z31.SinX*z32.X2SinX*z33.CosX)...
                    +(T23/T33)*(z31.CosX*z32.CosX*z33.X2CosX-z31.CosX*z32.SinX*z33.X2SinX-z31.SinX*z32.CosX*z33.X2SinX-z31.SinX*z32.SinX*z33.X2CosX)...
                    +((T21/T31)+(T22/T32))*(z31.XCosX*z32.XCosX*z33.CosX-z31.XCosX*z32.XSinX*z33.SinX-z31.XSinX*z32.XCosX*z33.SinX-z31.XSinX*z32.XSinX*z33.CosX)...
                    +((T22/T32)+(T23/T33))*(z31.CosX*z32.XCosX*z33.XCosX-z31.CosX*z32.XSinX*z33.XSinX-z31.SinX*z32.XCosX*z33.XSinX-z31.SinX*z32.XSinX*z33.XCosX)...
                    +((T21/T31)+(T23/T33))*(z31.XCosX*z32.CosX*z33.XCosX-z31.XCosX*z32.SinX*z33.XSinX-z31.XSinX*z32.CosX*z33.XSinX-z31.XSinX*z32.SinX*z33.XCosX);
            end
        end
        function e=YThS(obj)
            % calc E[(y*theta*sin(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.Mu(2)...
                    *((obj.Mu(3)*sin(obj.Mu(3))+obj.Sigma2(3,3)*cos(obj.Mu(3))) * exp(-0.5*obj.Sigma2(3,3)));
            else

                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T21/T31)*(-z31.X2SinX*z32.SinX*z33.SinX+z31.X2SinX*z32.CosX*z33.CosX+z31.X2CosX*z32.SinX*z33.CosX+z31.X2CosX*z32.CosX*z33.SinX)...
                    +(T22/T32)*(-z31.SinX*z32.X2SinX*z33.SinX+z31.SinX*z32.X2CosX*z33.CosX+z31.CosX*z32.X2SinX*z33.CosX+z31.CosX*z32.X2CosX*z33.SinX)...
                    +(T23/T33)*(-z31.SinX*z32.SinX*z33.X2SinX+z31.SinX*z32.CosX*z33.X2CosX+z31.CosX*z32.SinX*z33.X2CosX+z31.CosX*z32.CosX*z33.X2SinX)...
                    +((T21/T31)+(T22/T32))*(-z31.XSinX*z32.XSinX*z33.SinX+z31.XSinX*z32.XCosX*z33.CosX+z31.XCosX*z32.XSinX*z33.CosX+z31.XCosX*z32.XCosX*z33.SinX)...
                    +((T22/T32)+(T23/T33))*(-z31.SinX*z32.XSinX*z33.XSinX+z31.SinX*z32.XCosX*z33.XCosX+z31.CosX*z32.XSinX*z33.XCosX+z31.CosX*z32.XCosX*z33.XSinX)...
                    +((T21/T31)+(T23/T33))*(-z31.XSinX*z32.SinX*z33.XSinX+z31.XSinX*z32.CosX*z33.XCosX+z31.XCosX*z32.SinX*z33.XCosX+z31.XCosX*z32.CosX*z33.XSinX);
            end
        end


        function [T11,T12,T13,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj)
            T11 = obj.T(1,1);
            T12 = obj.T(1,2);
            T13 = obj.T(1,3);
            T21 = obj.T(2,1);
            T22 = obj.T(2,2);
            T23 = obj.T(2,3);
            T31 = obj.T(3,1);
            T32 = obj.T(3,2);
            T33 = obj.T(3,3);
            z31 = obj.Z31;
            z32 = obj.Z32;
            z33 = obj.Z33;
        end
    end
end
