classdef InEKF < handle   
    properties
        mu;                 % Pose Mean(3,3)
        Sigma;              % Pose Sigma(3,3)
        gfun;               % Motion model function
        mu_pred;             % Mean after prediction step
        Sigma_pred;          % Sigma after prediction step
        mu_cart;
        sigma_cart;
        M;                  % Motion model noise
        Q;                  % Sensor Noise
        H;
        A;
    end
    
    methods
        function obj = InEKF(sys, init)
            obj.gfun = sys.gfun;
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            
            obj.M = sys.M;
            obj.Q = [10 0 ; 0 10];
            obj.A = eye(3);
%             obj.Q = diag([0.015^2, 0.01^2, 0.01^2]);
%             obj.N = diag([0.5^2; 0.5^2]);
            
            
        end
        
        function prediction(obj, u)
            %Formulate Adjoint function to be used in propagation


            AdjX = [obj.mu(1:2,1:2), [obj.mu(2,3); -obj.mu(1,3)]; 0 0 1]; %

            %Convert motion command into lie algebra element to pass in to
            %propagation
            MU = [ obj.mu(1,3), obj.mu(2,3), atan2(obj.mu(2,1),obj.mu(1,1)) ];
            X = obj.posemat(obj.gfun( MU , u));
            U = logm(obj.mu \ X);
            obj.propagation(U , AdjX);
       
        end
        
        function propagation(obj, u, AdjX)
            % SE(2) propagation model; the input is u \in se(2) plus noise
            % propagate mean
            obj.mu_pred = obj.mu * expm(u);
        
            % propagate covariance
            obj.Sigma_pred = obj.A * obj.Sigma * obj.A' + AdjX * obj.M(u) * AdjX';
        end
        
        function correction(obj, Y, Y2, landmark_ids)
            global FIELDINFO;        
            landmark_x = FIELDINFO.MARKER_X_POS(landmark_ids(1));
            landmark_y = FIELDINFO.MARKER_Y_POS(landmark_ids(1));       
            landmark_x2 = FIELDINFO.MARKER_X_POS(landmark_ids(2));
            landmark_y2 = FIELDINFO.MARKER_Y_POS(landmark_ids(2));
            Y = Y';
            Y2 = Y2';
            
%             H1 = [  landmark_y, -1,  0;
%                    -landmark_x,  0, -1 ];
%             H2 = [  landmark_y2, -1,  0;
%                    -landmark_x2,  0, -1 ];
%                
            H1 = [-1,  0, landmark_y;
                   0, -1, -landmark_x];
            H2 = [-1,  0, landmark_y2;
                   0, -1, -landmark_x2];
            H = [H1;H2];
            
            
            N = obj.mu_pred * blkdiag(obj.Q,0) * obj.mu_pred'; N = N(1:2,1:2);
            
            S1 = H1 * obj.Sigma_pred * H1' + N; %3*3
            S2 = H2 * obj.Sigma_pred * H2' + N; %3*3
            Kgain1 = obj.Sigma * H1' / S1;
            Kgain2 = obj.Sigma * H2' / S2;
            Kgain = [Kgain1 Kgain2];
            b1 = [landmark_x, landmark_y, 1]';
            b2 = [landmark_x2, landmark_y2, 1]';
            
            nu1 = (obj.mu_pred * Y -  b1); nu1 = nu1(1:2); %%%!!!!!
            nu2 = (obj.mu_pred * Y2 - b2); nu2 = nu2(1:2);
            
            delta1 = obj.wedge(Kgain1 * nu1);
            delta2 = obj.wedge(Kgain2 * nu2);
            delta = delta1 + delta2;
            
            obj.mu = expm(delta) * obj.mu_pred;
            obj.Sigma = (eye(3) - Kgain * H) * obj.Sigma_pred * (eye(3) - Kgain * H)'...
                            + Kgain * blkdiag(N,N) * Kgain';
        end
        
        function H = posemat(obj,state)
            x = state(1);
            y = state(2);
            h = state(3);
            % construct a SE(2) matrix element
            H = [...
                cos(h) -sin(h) x;
                sin(h)  cos(h) y;
                     0       0 1];
        end
        
        function xhat = wedge(obj,x)
            % wedge operation for se(2) to put an R^3 vector into the Lie
            % algebra basis.
            G1 = [0    -1     0
                1     0     0
                0     0     0] ; % omega
            G2 = [0     0     1
                0     0     0
                0     0     0] ; % v_1
            G3 = [0     0     0
                0     0     1
                0     0     0] ; % v_2
            xhat = G1 * x(1) + G2 * x(2) + G3 * x(3);
        end
    end
end
