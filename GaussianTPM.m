% 
% compute E[f(x)],x~N(mu sigma2)
% input: mu,sigma (scalar) 
% include Mixed-trigonometric-Polynominal Moments
% E[\x^a1 cos^a2(x) sin^a3(x)],
% refference:https://arxiv.org/pdf/2101.12490.pdf
classdef GaussianTPM
    properties
        mu
        sigma2
    end
    methods
        function obj = GaussianTPM(mu,sigma2)
            if nargin == 2
                obj.mu = mu;
                obj.sigma2 = sigma2;
            end
        end
        function e = X(obj)
            e = obj.mu;
        end
        function e = X2(obj)
            % E[x^2]
            e = obj.sigma2+obj.mu^2;
        end
        function e = CosX(obj)
            % E[cosx]
            e = cos(obj.mu)*exp(-0.5*obj.sigma2);
        end
        function e = SinX(obj)
            % E[sinx]
            e = sin(obj.mu)*exp(-0.5*obj.sigma2);
        end
        function e = XCosX(obj)
            % E[xcosx]
            e = (obj.mu*cos(obj.mu)-obj.sigma2*sin(obj.mu)) * exp(-0.5*obj.sigma2);
        end
        function e = XSinX(obj)
            % E[xsinx]
            e = (obj.mu*sin(obj.mu)+obj.sigma2*cos(obj.mu)) * exp(-0.5*obj.sigma2);
        end 
        function e = X2CosX(obj)
            % E[x^2 cosx]
            e = ((obj.sigma2+obj.mu^2-obj.sigma2^2)*cos(obj.mu)-2*obj.mu*obj.sigma2*sin(obj.mu))*exp(-0.5*obj.sigma2);
        end 
        function e = X2SinX(obj)
            % E[x^2 sinx]
            e =((obj.sigma2+obj.mu^2-obj.sigma2^2)*sin(obj.mu)+2*obj.mu*obj.sigma2*cos(obj.mu))*exp(-0.5*obj.sigma2);
        end 
    end
end
