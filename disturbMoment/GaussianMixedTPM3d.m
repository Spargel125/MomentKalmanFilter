classdef GaussianMixedTPM3d < GaussianTPM
    properties
        Mu,Sigma2
        x,y,th
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
                obj.x = GaussianTPM(obj.Mu(1),obj.Sigma2(1,1));
                obj.y = GaussianTPM(obj.Mu(2),obj.Sigma2(2,2));
                obj.th = GaussianTPM(obj.Mu(3),obj.Sigma2(3,3));
            end
        end
        function e = X(obj)
            e = obj.x.X;
        end
        function e = Y(obj)
            e = obj.y.X;
        end
        function e = Th(obj)
            e = obj.th.X;
        end
        function e = C(obj)
            e = obj.th.C;
        end
        function e = S(obj)
            e = obj.th.S;
        end
        function e = CS(obj)
            %             e = 0.5*sin(2*obj.Mu(3))*exp(-0.5*4*obj.Sigma2(3,3));
            e = obj.th.CS;
        end
        function e = SC(obj)
            e = obj.th.CS;
        end
        function e = CC(obj)
            %             e = 0.5*(cos(2*obj.Mu(3))*exp(-0.5*4*obj.Sigma2(3,3))+1);
            e = obj.th.CC;
        end
        function e = SS(obj)
            %             e = 0.5*(1-cos(2*obj.Mu(3))*exp(-0.5*4*obj.Sigma2(3,3)));
            e = obj.th.SS;
        end

        function e = XX(obj)
            % E[x^2]
            %             e = obj.Sigma2(1,1)+obj.Mu(1)^2;
            e = obj.x.XX;
        end
        function e = YY(obj)
            % E[y^2]
            %             e = obj.Sigma2(2,2)+obj.Mu(2)^2;
            e = obj.y.XX;
        end
        function e = ThTh(obj)
            % E[y^2]
            %             e = obj.Sigma2(3,3)+obj.Mu(3)^2;
            e = obj.th.XX;
        end

        function e = ThC(obj)
            %             e = (obj.Mu(3)*cos(obj.Mu(3))-obj.Sigma2(3,3)*sin(obj.Mu(3))) * exp(-0.5*obj.Sigma2(3,3));
            e = obj.th.XC;

        end
        function e = ThS(obj)
            %             e = (obj.Mu(3)*sin(obj.Mu(3))+obj.Sigma2(3,3)*cos(obj.Mu(3))) * exp(-0.5*obj.Sigma2(3,3));
            e = obj.th.XS;
        end

        function e = XY(obj)
            if obj.Sigma2(1,2)==0
                e = obj.x.X*obj.y.X;
            else
                [T11,T12,T13,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*((T21/T31)*z31.XX+(T22/T32)*z31.X*z32.X+(T23/T33)*z31.X*z33.X)...
                    +(T12/T32)*((T21/T31)*z31.X*z32.X+(T22/T32)*z32.XX+(T23/T33)*z32.X*z33.X)...
                    +(T13/T33)*((T21/T31)*z33.X*z31.X+(T22/T32)*z33.X*z32.X+(T23/T33)*z33.XX);
            end
        end
        function e = XTh(obj)
            if obj.Sigma2(1,3)==0
                e = obj.x.X*obj.th.X;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(z31.XX+z31.X*z32.X+z31.X*z33.X)...
                    +(T12/T32)*(z31.X*z32.X+z32.XX+z32.X*z33.X)...
                    +(T13/T33)*(z33.X*z31.X+z33.X*z32.X+z33.XX);
            end
        end
        function e = YTh(obj)
            if obj.Sigma2(2,3)==0
                e = obj.y.X*obj.th.X;
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T21/T31)*(z31.XX+z31.X*z32.X+z31.X*z33.X)...
                    +(T22/T32)*(z31.X*z32.X+z32.XX+z32.X*z33.X)...
                    +(T23/T33)*(z33.X*z31.X+z33.X*z32.X+z33.XX);
            end
        end

        function e = XC(obj)
            % calc E[(x*cos(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.C;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)*(z31.XC*z32.C*z33.C-z31.XC*z32.S*z33.S-z31.XS*z32.C*z33.S-z31.XS*z32.S*z33.C)...
                    +(T12/T32)*(z31.C*z32.XC*z33.C-z31.C*z32.XS*z33.S-z31.S*z32.XC*z33.S-z31.S*z32.XS*z33.C)...
                    +(T13/T33)*(z31.C*z32.C*z33.XC-z31.C*z32.S*z33.XS-z31.S*z32.C*z33.XS-z31.S*z32.S*z33.XC);
            end
        end
        function e = XS(obj)
            % calc E[(x*sin(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.S;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)*(-z31.XS*z32.S*z33.S+z31.XS*z32.C*z33.C+z31.XC*z32.S*z33.C+z31.XC*z32.C*z33.S)...
                    +(T12/T32)*(-z31.S*z32.XS*z33.S+z31.S*z32.XC*z33.C+z31.C*z32.XS*z33.C+z31.C*z32.XC*z33.S)...
                    +(T13/T33)*(-z31.S*z32.S*z33.XS+z31.S*z32.C*z33.XC+z31.C*z32.S*z33.XC+z31.C*z32.C*z33.XS);
            end
        end
        function e = YC(obj)
            % calc E[(y*cos(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.C;
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)*(z31.XC*z32.C*z33.C-z31.XC*z32.S*z33.S-z31.XS*z32.C*z33.S-z31.XS*z32.S*z33.C)...
                    +(T22/T32)*(z31.C*z32.XC*z33.C-z31.C*z32.XS*z33.S-z31.S*z32.XC*z33.S-z31.S*z32.XS*z33.C)...
                    +(T23/T33)*(z31.C*z32.C*z33.XC-z31.C*z32.S*z33.XS-z31.S*z32.C*z33.XS-z31.S*z32.S*z33.XC);
            end
        end
        function e = YS(obj)
            % calc E[(y*sin(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.S;
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)*(-z31.XS*z32.S*z33.S+z31.XS*z32.C*z33.C+z31.XC*z32.S*z33.C+z31.XC*z32.C*z33.S)...
                    +(T22/T32)*(-z31.S*z32.XS*z33.S+z31.S*z32.XC*z33.C+z31.C*z32.XS*z33.C+z31.C*z32.XC*z33.S)...
                    +(T23/T33)*(-z31.S*z32.S*z33.XS+z31.S*z32.C*z33.XC+z31.C*z32.S*z33.XC+z31.C*z32.C*z33.XS);
            end
        end
        function e = XXC(obj)
            % calc E[(x*x*cos(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.XX*obj.th.C;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)^2*(z31.XXC*z32.C*z33.C-z31.XXC*z32.S*z33.S-z31.XXS*z32.C*z33.S-z31.XXS*z32.S*z33.C)...
                    +(T12/T32)^2*(z31.C*z32.XXC*z33.C-z31.C*z32.XXS*z33.S-z31.S*z32.XXC*z33.S-z31.S*z32.XXS*z33.C)...
                    +(T13/T33)^2*(z31.C*z32.C*z33.XXC-z31.C*z32.S*z33.XXS-z31.S*z32.C*z33.XXS-z31.S*z32.S*z33.XXC)...
                    +2*(T11/T31)*(T12/T32)*(z31.XC*z32.XC*z33.C-z31.XC*z32.XS*z33.S-z31.XS*z32.XC*z33.S-z31.XS*z32.XS*z33.C)...
                    +2*(T12/T32)*(T13/T33)*(z31.C*z32.XC*z33.XC-z31.C*z32.XS*z33.XS-z31.S*z32.XC*z33.XS-z31.S*z32.XS*z33.XC)...
                    +2*(T13/T33)*(T11/T31)*(z31.XC*z32.C*z33.XC-z31.XC*z32.S*z33.XS-z31.XS*z32.C*z33.XS-z31.XS*z32.S*z33.XC);
            end
        end
        function e = XCC(obj)
            % calc E[(x**cos^2(theta))]=E[x**(cos(2theta)+1)/2]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.CC;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/4)*(tmp.X+tmp.XC);
            end
        end
        function e = XSS(obj)
            % calc E[(x*x*sin^2(theta))]=E[x*(1-cos(2theta))/2]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.SS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/4)*(tmp.X-tmp.XC);
            end
        end
        function e = XCS(obj)
            % calc E[(x*cos(theta)*sin(theta)]=E[x*sin(2theta)/2]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.CS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/4)*(tmp.XS);
            end
        end
        function e = XSC(obj)
            e = obj.XCS;
        end
        function e = XXCS(obj)
            % calc E[(xx*cos(theta)*sin(theta)]=E[xx*sin(2theta)/2]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.XX*obj.th.CS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.XXS);
            end
        end
        function e = XXSC(obj)
            e = obj.XXCS;
        end

        function e = YCC(obj)
            % calc E[(y*cos^2(theta))]=E[y*(cos(2theta)+1)/2]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.CC;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/4)*(tmp.Y+tmp.YC);
            end
        end
        function e = YSS(obj)
            % calc E[(y*sin^2(theta))]=E[y*(1-cos(2theta))/2]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.SS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/4)*(tmp.Y-tmp.YC);
            end
        end
        function e = YCS(obj)
            % calc E[(y*cos(theta)*sin(theta)]=E[y*sin(2theta)/2]
            if obj.Sigma2(2,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.CS;
            else
                %y=2y,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/4)*(tmp.YS);
            end
        end
        function e = YSC(obj)
            e = obj.YCS;
        end
        function e = YYCS(obj)
            % calc E[(yy*cos(theta)*sin(theta)]=E[yy*sin(2theta)/2]
            if obj.Sigma2(2,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.XX*obj.th.CS;
            else
                %y=2y,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.YYS);
            end
        end
        function e = YYSC(obj)
            e = obj.YYCS;
        end

        function e = XYCS(obj)
            % calc E[(x*y*cos(theta)*sin(theta)]=E[x*y*sin(2theta)/2]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = obj.x.X*obj.y.X*obj.th.CS;
            else
                %x=2X,y=2y,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.XYS);
            end
        end
        function e = XYSC(obj)
            e = obj.XYCS;
        end
        function e = XYCC(obj)
            % calc E[(xy*cos^2(theta))]=E[xy*(cos(2theta)+1)/2]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = obj.x.X*obj.y.X*obj.th.CC;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.XY+tmp.XYC);
            end
        end
        function e = XYSS(obj)
            % calc E[(xy*sin^2(theta))]=E[xy*(1-cos(2theta))/2]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = obj.x.X*obj.y.X*obj.th.SS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.XY-tmp.XYC);
            end
        end

        function e = XXCC(obj)
            % calc E[(x*x*cos^2(theta))]=E[x*x*(cos(2theta)+1)/2]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.XX*obj.th.CC;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.XX+tmp.XXC);
            end
        end
        function e = XXSS(obj)
            % calc E[(x*x*sin^2(theta))]=E[x*x*(1-cos(2theta))/2]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.XX*obj.th.SS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.XX-tmp.XXC);
            end
        end
        function e = YYCC(obj)
            % calc E[(y*y*cos^2(theta))]=E[y*y*(cos(2theta)+1)/2]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.XX*obj.th.CC;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.YY+tmp.YYC);
            end
        end
        function e = YYSS(obj)
            % calc E[(y*y*sin^2(theta))]=E[y*y*(1-cos(2theta))/2]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.XX*obj.th.SS;
            else
                %x=2x,theta=2thetaと置きなおして再帰的に呼び出す
                tmp = GaussianMixedTPM3d(2*obj.Mu,4*obj.Sigma2);
                e = (1/8)*(tmp.YY-tmp.YYC);
            end
        end

        function e = XXS(obj)
            % calc E[(x*x*sin(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.XX*obj.th.S;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T11/T31)^2*(-z31.XXS*z32.S*z33.S+z31.XXS*z32.C*z33.C+z31.XXC*z32.S*z33.C+z31.XXC*z32.C*z33.S)...
                    +(T12/T32)^2*(-z31.S*z32.XXS*z33.S+z31.S*z32.XXC*z33.C+z31.C*z32.XXS*z33.C+z31.C*z32.XXC*z33.S)...
                    +(T13/T33)^2*(-z31.S*z32.S*z33.XXS+z31.S*z32.C*z33.XXC+z31.C*z32.S*z33.XXC+z31.C*z32.C*z33.XXS)...
                    +2*(T11/T31)*(T12/T32)*(-z31.XS*z32.XS*z33.S+z31.XS*z32.XC*z33.C+z31.XC*z32.XS*z33.C+z31.XC*z32.XC*z33.S)...
                    +2*(T12/T32)*(T13/T33)*(-z31.S*z32.XS*z33.XS+z31.S*z32.XC*z33.XC+z31.C*z32.XS*z33.XC+z31.C*z32.XC*z33.XS)...
                    +2*(T13/T33)*(T11/T31)*(-z31.XS*z32.S*z33.XS+z31.XS*z32.C*z33.XC+z31.XC*z32.S*z33.XC+z31.XC*z32.C*z33.XS);
            end
        end
        function e = YYC(obj)
            % calc E[(y*y*cos(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.XX*obj.th.C;
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)^2*(z31.XXC*z32.C*z33.C-z31.XXC*z32.S*z33.S-z31.XXS*z32.C*z33.S-z31.XXS*z32.S*z33.C)...
                    +(T22/T32)^2*(z31.C*z32.XXC*z33.C-z31.C*z32.XXS*z33.S-z31.S*z32.XXC*z33.S-z31.S*z32.XXS*z33.C)...
                    +(T23/T33)^2*(z31.C*z32.C*z33.XXC-z31.C*z32.S*z33.XXS-z31.S*z32.C*z33.XXS-z31.S*z32.S*z33.XXC)...
                    +2*(T21/T31)*(T22/T32)*(z31.XC*z32.XC*z33.C-z31.XC*z32.XS*z33.S-z31.XS*z32.XC*z33.S-z31.XS*z32.XS*z33.C)...
                    +2*(T22/T32)*(T23/T33)*(z31.C*z32.XC*z33.XC-z31.C*z32.XS*z33.XS-z31.S*z32.XC*z33.XS-z31.S*z32.XS*z33.XC)...
                    +2*(T23/T33)*(T21/T31)*(z31.XC*z32.C*z33.XC-z31.XC*z32.S*z33.XS-z31.XS*z32.C*z33.XS-z31.XS*z32.S*z33.XC);
            end
        end
        function e = YYS(obj)
            % calc E[(y*y*sin(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.XX*obj.th.S;
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e =  (T21/T31)^2*(-z31.XXS*z32.S*z33.S+z31.XXS*z32.C*z33.C+z31.XXC*z32.S*z33.C+z31.XXC*z32.C*z33.S)...
                    +(T22/T32)^2*(-z31.S*z32.XXS*z33.S+z31.S*z32.XXC*z33.C+z31.C*z32.XXS*z33.C+z31.C*z32.XXC*z33.S)...
                    +(T23/T33)^2*(-z31.S*z32.S*z33.XXS+z31.S*z32.C*z33.XXC+z31.C*z32.S*z33.XXC+z31.C*z32.C*z33.XXS)...
                    +2*(T21/T31)*(T22/T32)*(-z31.XS*z32.XS*z33.S+z31.XS*z32.XC*z33.C+z31.XC*z32.XS*z33.C+z31.XC*z32.XC*z33.S)...
                    +2*(T22/T32)*(T23/T33)*(-z31.S*z32.XS*z33.XS+z31.S*z32.XC*z33.XC+z31.C*z32.XS*z33.XC+z31.C*z32.XC*z33.XS)...
                    +2*(T23/T33)*(T21/T31)*(-z31.XS*z32.S*z33.XS+z31.XS*z32.C*z33.XC+z31.XC*z32.S*z33.XC+z31.XC*z32.C*z33.XS);
            end
        end
        function e = XYC(obj)
            % calc E[(x*y*cos(theta))]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = obj.x.X*obj.y.X*obj.th.C;
            else
                [T11,T12,T13,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(T21/T31)*(z31.XXC*z32.C*z33.C-z31.XXC*z32.S*z33.S-z31.XXS*z32.C*z33.S-z31.XXS*z32.S*z33.C)...
                    +(T12/T32)*(T22/T32)*(z31.C*z32.XXC*z33.C-z31.C*z32.XXS*z33.S-z31.S*z32.XXC*z33.S-z31.S*z32.XXS*z33.C)...
                    +(T13/T33)*(T23/T33)*(z31.C*z32.C*z33.XXC-z31.C*z32.S*z33.XXS-z31.S*z32.C*z33.XXS-z31.S*z32.S*z33.XXC)...
                    +((T11/T31)*(T22/T32)+(T12/T32)*(T21/T31))*(z31.XC*z32.XC*z33.C-z31.XC*z32.XS*z33.S-z31.XS*z32.XC*z33.S-z31.XS*z32.XS*z33.C)...
                    +((T12/T32)*(T23/T33)+(T13/T33)*(T22/T32))*(z31.C*z32.XC*z33.XC-z31.C*z32.XS*z33.XS-z31.S*z32.XC*z33.XS-z31.S*z32.XS*z33.XC)...
                    +((T13/T33)*(T21/T31)+(T11/T31)*(T23/T33))*(z31.XC*z32.C*z33.XC-z31.XC*z32.S*z33.XS-z31.XS*z32.C*z33.XS-z31.XS*z32.S*z33.XC);
            end
        end
        function e = XYS(obj)
            % calc E[(x*y*sin(theta))]
            if isdiag(obj.Sigma2) == true %Sigma2が対角→x,y,thetaが独立
                e = obj.x.X*obj.y.X*obj.th.S;
            else
                [T11,T12,T13,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(T21/T31)*(-z31.XXS*z32.S*z33.S+z31.XXS*z32.C*z33.C+z31.XXC*z32.S*z33.C+z31.XXC*z32.C*z33.S)...
                    +(T12/T32)*(T22/T32)*(-z31.S*z32.XXS*z33.S+z31.S*z32.XXC*z33.C+z31.C*z32.XXS*z33.C+z31.C*z32.XXC*z33.S)...
                    +(T13/T33)*(T23/T33)*(-z31.S*z32.S*z33.XXS+z31.S*z32.C*z33.XXC+z31.C*z32.S*z33.XXC+z31.C*z32.C*z33.XXS)...
                    +((T11/T31)*(T22/T32)+(T12/T32)*(T21/T31))*(-z31.XS*z32.XS*z33.S+z31.XS*z32.XC*z33.C+z31.XC*z32.XS*z33.C+z31.XC*z32.XC*z33.S)...
                    +((T12/T32)*(T23/T33)+(T13/T33)*(T22/T32))*(-z31.S*z32.XS*z33.XS+z31.S*z32.XC*z33.XC+z31.C*z32.XS*z33.XC+z31.C*z32.XC*z33.XS)...
                    +((T13/T33)*(T21/T31)+(T11/T31)*(T23/T33))*(-z31.XS*z32.S*z33.XS+z31.XS*z32.C*z33.XC+z31.XC*z32.S*z33.XC+z31.XC*z32.C*z33.XS);
            end
        end
        function e=XThC(obj)
            % calc E[(x*theta*cos(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.XC;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(z31.XXC*z32.C*z33.C-z31.XXC*z32.S*z33.S-z31.XXS*z32.C*z33.S-z31.XXS*z32.S*z33.C)...
                    +(T12/T32)*(z31.C*z32.XXC*z33.C-z31.C*z32.XXS*z33.S-z31.S*z32.XXC*z33.S-z31.S*z32.XXS*z33.C)...
                    +(T13/T33)*(z31.C*z32.C*z33.XXC-z31.C*z32.S*z33.XXS-z31.S*z32.C*z33.XXS-z31.S*z32.S*z33.XXC)...
                    +((T11/T31)+(T12/T32))*(z31.XC*z32.XC*z33.C-z31.XC*z32.XS*z33.S-z31.XS*z32.XC*z33.S-z31.XS*z32.XS*z33.C)...
                    +((T12/T32)+(T13/T33))*(z31.C*z32.XC*z33.XC-z31.C*z32.XS*z33.XS-z31.S*z32.XC*z33.XS-z31.S*z32.XS*z33.XC)...
                    +((T11/T31)+(T13/T33))*(z31.XC*z32.C*z33.XC-z31.XC*z32.S*z33.XS-z31.XS*z32.C*z33.XS-z31.XS*z32.S*z33.XC);
            end
        end
        function e=XThS(obj)
            % calc E[(x*theta*sin(theta))]
            if obj.Sigma2(1,3) == 0 % xとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.x.X*obj.th.XS;
            else
                [T11,T12,T13,~,~,~,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T11/T31)*(-z31.XXS*z32.S*z33.S+z31.XXS*z32.C*z33.C+z31.XXC*z32.S*z33.C+z31.XXC*z32.C*z33.S)...
                    +(T12/T32)*(-z31.S*z32.XXS*z33.S+z31.S*z32.XXC*z33.C+z31.C*z32.XXS*z33.C+z31.C*z32.XXC*z33.S)...
                    +(T13/T33)*(-z31.S*z32.S*z33.XXS+z31.S*z32.C*z33.XXC+z31.C*z32.S*z33.XXC+z31.C*z32.C*z33.XXS)...
                    +((T11/T31)+(T12/T32))*(-z31.XS*z32.XS*z33.S+z31.XS*z32.XC*z33.C+z31.XC*z32.XS*z33.C+z31.XC*z32.XC*z33.S)...
                    +((T12/T32)+(T13/T33))*(-z31.S*z32.XS*z33.XS+z31.S*z32.XC*z33.XC+z31.C*z32.XS*z33.XC+z31.C*z32.XC*z33.XS)...
                    +((T11/T31)+(T13/T33))*(-z31.XS*z32.S*z33.XS+z31.XS*z32.C*z33.XC+z31.XC*z32.S*z33.XC+z31.XC*z32.C*z33.XS);
            end
        end
        function e=YThC(obj)
            % calc E[(y*theta*cos(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.XC;
            else

                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T21/T31)*(z31.XXC*z32.C*z33.C-z31.XXC*z32.S*z33.S-z31.XXS*z32.C*z33.S-z31.XXS*z32.S*z33.C)...
                    +(T22/T32)*(z31.C*z32.XXC*z33.C-z31.C*z32.XXS*z33.S-z31.S*z32.XXC*z33.S-z31.S*z32.XXS*z33.C)...
                    +(T23/T33)*(z31.C*z32.C*z33.XXC-z31.C*z32.S*z33.XXS-z31.S*z32.C*z33.XXS-z31.S*z32.S*z33.XXC)...
                    +((T21/T31)+(T22/T32))*(z31.XC*z32.XC*z33.C-z31.XC*z32.XS*z33.S-z31.XS*z32.XC*z33.S-z31.XS*z32.XS*z33.C)...
                    +((T22/T32)+(T23/T33))*(z31.C*z32.XC*z33.XC-z31.C*z32.XS*z33.XS-z31.S*z32.XC*z33.XS-z31.S*z32.XS*z33.XC)...
                    +((T21/T31)+(T23/T33))*(z31.XC*z32.C*z33.XC-z31.XC*z32.S*z33.XS-z31.XS*z32.C*z33.XS-z31.XS*z32.S*z33.XC);
            end
        end
        function e=YThS(obj)
            % calc E[(y*theta*sin(theta))]
            if obj.Sigma2(2,3) == 0 % yとθの共分散が0→独立なのでそれぞれの確率モーメントの積で返す
                e = obj.y.X*obj.th.XS;
            else
                [~,~,~,T21,T22,T23,T31,T32,T33,z31,z32,z33] = variableExpansion(obj);
                e = (T21/T31)*(-z31.XXS*z32.S*z33.S+z31.XXS*z32.C*z33.C+z31.XXC*z32.S*z33.C+z31.XXC*z32.C*z33.S)...
                    +(T22/T32)*(-z31.S*z32.XXS*z33.S+z31.S*z32.XXC*z33.C+z31.C*z32.XXS*z33.C+z31.C*z32.XXC*z33.S)...
                    +(T23/T33)*(-z31.S*z32.S*z33.XXS+z31.S*z32.C*z33.XXC+z31.C*z32.S*z33.XXC+z31.C*z32.C*z33.XXS)...
                    +((T21/T31)+(T22/T32))*(-z31.XS*z32.XS*z33.S+z31.XS*z32.XC*z33.C+z31.XC*z32.XS*z33.C+z31.XC*z32.XC*z33.S)...
                    +((T22/T32)+(T23/T33))*(-z31.S*z32.XS*z33.XS+z31.S*z32.XC*z33.XC+z31.C*z32.XS*z33.XC+z31.C*z32.XC*z33.XS)...
                    +((T21/T31)+(T23/T33))*(-z31.XS*z32.S*z33.XS+z31.XS*z32.C*z33.XC+z31.XC*z32.S*z33.XC+z31.XC*z32.C*z33.XS);
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
