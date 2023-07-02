classdef GaussianMixedTPM2d < GaussianTPM
    properties
        Mu,Sigma2
        T,Lambda
        Y1,Y2
        L1,L2
        M1,M2
    end
    methods
        function obj = GaussianMixedTPM2d(Mu,Sigma2)
            if nargin == 2
                if size(Mu,1) ~= 2
                    error('Not 2-dimension, check dimension of mu or sigma')
                end
                obj.Mu = Mu;
                obj.Sigma2 = Sigma2;
                [obj.T,obj.Lambda] = eig(Sigma2);
                T21 = obj.T(2,1);
                T22 = obj.T(2,2);
                mu_y = obj.T'*obj.Mu;
                mu_y1 = mu_y(1);
                mu_y2 = mu_y(2);
                lambda1 = obj.Lambda(1,1);
                lambda2 = obj.Lambda(2,2);
                % 確率分布のクラス定義
                obj.Y1 = GaussianTPM(mu_y1,lambda1);
                obj.Y2 = GaussianTPM(mu_y2,lambda2);
                obj.L1 = GaussianTPM(T21*mu_y1,T21^2*lambda1);
                obj.L2 = GaussianTPM(T22*mu_y2,T22^2*lambda2);
                obj.M1 = GaussianTPM(2*T21*mu_y1,(2*T21)^2*lambda1);
                obj.M2 = GaussianTPM(2*T22*mu_y2,(2*T22)^2*(lambda2));
            end
        end
        function e = XTh(obj)
            % E[x*theta]
            if obj.Sigma2(1,2)==0 %２変数が独立の場合，それぞれの積で計算する．
                e = obj.Mu(1)...
                    *obj.Mu(2);
            else
                [T11,T12,T21,T22,y1,y2,~] = variableExpansion(obj);
                e = T11*T21*y1.X2+(T11*T22+T12*T21)*y1.X*y2.X+T12*T22*y2.X2;

            end
        end
        function e = XC(obj)
            if obj.Sigma2(1,2)==0 %２変数が独立の場合，それぞれの積で計算する．
                e = obj.Mu(1)...
                    *(cos(obj.Mu(2))*exp(-0.5*obj.Sigma2(2,2)));
            else
            % E[x*cos(theta)]
            [T11,T12,T21,T22,~,~,l1,l2,~] = variableExpansion(obj);
            e = T11/T21*(l1.XCosX*l2.CosX-l1.XSinX*l2.SinX)+T12/T22*(l1.CosX*l2.XCosX-l1.SinX*l2.XSinX);
            end
        end
        function e = XS(obj)
            % E[x*sin(theta)]
            if obj.Sigma2(1,2)==0 %２変数が独立の場合，それぞれの積で計算する．
                e = obj.Mu(1)...
                    *(sin(obj.Mu(2))*exp(-0.5*obj.Sigma2(2,2)));
            else
            [T11,T12,T21,T22,~,~,l1,l2,~] = variableExpansion(obj);
            e = T11/T21*(l1.XSinX*l2.CosX+l1.XCosX*l2.SinX)+T12/T22*(l1.SinX*l2.XCosX+l1.CosX*l2.XSinX);
            end
        end
        function e = XCS(obj)
            % E[x*cos(theta)*sin(theta)]
            if obj.Sigma2(1,2)==0 %２変数が独立の場合，それぞれの積で計算する．
                e = obj.Mu(1)...
                    *(0.5*sin(2*obj.Mu(2))*exp(-0.5*4*obj.Sigma2(2,2)));
            else
            [T11,T12,T21,T22,~,~,~,~,m1,m2] = variableExpansion(obj);
            e = T11/(4*T21)*(m1.XSinX*m2.CosX+m1.XCosX*m2.SinX)+T12/(4*T22)*(m1.SinX*m2.XCosX+m1.CosX*m2.XSinX);
            end
        end

        function [T11,T12,T21,T22,y1,y2,l1,l2,m1,m2] = variableExpansion(obj)
            T11 = obj.T(1,1);
            T12 = obj.T(1,2);
            T21 = obj.T(2,1);
            T22 = obj.T(2,2);
            y1 = obj.Y1;
            y2 = obj.Y2;
            m1 = obj.M1;
            m2 = obj.M2;
            l1 = obj.L1;
            l2 = obj.L2;
        end
    end
end
