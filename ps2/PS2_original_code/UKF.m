classdef UKF < handle
    properties
        mu;             % Pose Mean
        Sigma;          % Pose Covariance
        gfun;           % Motion Model Function
        hfun;           % Measruement Model Function
        M;              % Motion model noise(dynamical and function of input)
        Q;              % Sensor Noise
        kappa_g;
        mu_pred;
        Sigma_pred;
        n;
        X;
        w;
        Y;
    end
    
    methods
        function obj = UKF(sys, init)
            % motion model
            obj.gfun = sys.gfun;
            % measurement model
            obj.hfun = sys.hfun;
            % motion noise covariance
            obj.M = sys.M;
            % measurement noise covariance
            obj.Q = sys.Q;
            obj.kappa_g = init.kappa_g;
            % initial mean and covariance
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            
       
        end
        
        function prediction(obj, u)
            %Create Sigma Points
            sigma_point(obj, obj.mu, obj.Sigma, obj.kappa_g);       %Create W,X,Y,n,L
            obj.mu_pred = 0;
            obj.Sigma_pred = 0;
            % (Predicted) compute sample mean and covariance
            for j = 1 : 2 * obj.n + 1
                obj.Y(:,j) = obj.gfun(obj.X(:,j) , u);
                obj.mu_pred = obj.mu_pred + obj.w(j) * obj.Y(:,j);       %mu_pred
                diff = obj.Y(:,j) - obj.mu_pred;
%                 diff(1,:) = wrapToPi(diff(1,:));
                obj.Sigma_pred = obj.Sigma_pred + obj.w(j) * (diff) * (diff)';
            end
            obj.Sigma_pred = obj.Sigma_pred + obj.M(u); %sigma_pred
            
            %Update Covariance           

        end
        
        function correction(obj, z)
            global FIELDINFO;
            landmark_x = FIELDINFO.MARKER_X_POS(z(3));
            landmark_y = FIELDINFO.MARKER_Y_POS(z(3));
            
            sigma_point(obj , obj.mu_pred, obj.Sigma_pred, obj.kappa_g);
            %Calculate Innovation
            zhat = 0;
            S = 0;
            Cov_xz = 0;
            for j = 1 : 2 * obj.n + 1
                zhat = zhat + obj.w(j) * obj.hfun(landmark_x, landmark_y, obj.X(:,j)); % 2*1
            end
            z(1,:) = wrapToPi(z(1,:)); % 3*1
            nu = z(1:2) - zhat;
%             nu(2) = wrapToPi(nu(2));
            for j = 1 : 2 * obj.n + 1
                HX = obj.hfun(landmark_x, landmark_y, obj.X(:,j)); %2*1
                S = S + obj.w(j) * (HX - zhat) * (HX - zhat).';
                Cov_xz = Cov_xz + obj.w(j) * (obj.X(:,j) - obj.mu_pred) * (HX - zhat)';
            end           
            S = S + obj.Q; 
            %Calculate Kalman Gain
            Kgain = Cov_xz * inv(S);
            
            %Corrected!!!
            obj.mu = obj.mu_pred + Kgain * nu;
            obj.Sigma = obj.Sigma_pred - Kgain * S * Kgain';
%             obj.Sigma = (eye(3) - Kgain * S) * obj.Sigma_pred * (eye(3) - Kgain * S)' + Kgain * obj.Q * Kgain';
        end
        
        
        
        function sigma_point(obj, mean, cov, kappa)
            obj.n = numel(mean);
            L = sqrt(obj.n + kappa) * chol(cov,'lower');
            obj.Y = mean(:,ones(1, numel(mean)));
            obj.X = [mean,obj.Y + L, obj.Y - L]; %sigma points
            obj.w = zeros(2 * obj.n + 1,1);
            obj.w(1) = kappa / (obj.n + kappa);
            obj.w(2:end) = 0.5 / (obj.n + kappa);%calculate the weights
        end
    end
end