classdef PF < handle
    properties
        gfun;               %Motion model function
        hfun;               %Measurement Model Function       
        Q;                  %Sensor Noise
        M;                  %Motion Model Noise (dynamical and function of input)
        n;                  %Number of Particles
        particles;          %Pose of particle
        particle_weight;    %Particle Weight
        mu;
        Sigma;
    end
    
    methods
        function obj = PF(sys, init)
            % motion model
            obj.gfun = sys.gfun;
            % measurement model
            obj.hfun = sys.hfun;
            % motion noise covariance
            obj.M = sys.M;
            % measurement noise covariance
            obj.Q = sys.Q;
            % PF parameters
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            obj.n = init.n;
            obj.particles = init.particles;
            obj.particle_weight = init.particle_weight;
        end
        %% Prediction
        function prediction(obj, u)
            %Particles created with motion model function with input u
            %having noise(mu = u, noise cov = M(u))   
            for i = 1 : obj.n
                particles_noise = mvnrnd(u , obj.M(u));
                obj.particles(:,i) = obj.gfun(obj.particles(:,i), particles_noise );
                
            end
        end
        
        %% Correction
        function correction(obj, z)
            global FIELDINFO;
            landmark_x = FIELDINFO.MARKER_X_POS(z(3));
            landmark_y = FIELDINFO.MARKER_Y_POS(z(3));
            
            % update and normalize weights
            weight_plus = 0;
            for i = 1:obj.n
                 %Calculate measurement and difference in measurements
                 % We know here z(2) is an angle
                 % zhat is what we assume to be the correct mu, with input
                 % of mu_pred
                 
                 zhat = obj.hfun(landmark_x , landmark_y , obj.particles(:,i));
                 zhat(1) = wrapToPi(zhat(1));
                 %Use mvnpdf to get weight of corresponding measurement
                 w(i) = mvnpdf(z(1:2) , zhat , obj.Q);
                 before_normalized(i) = w(i) * obj.particle_weight(i); % weight = prev_weight * p(z|x_k)
                 weight_plus = weight_plus + before_normalized(i);
                 
            end
            obj.particle_weight = before_normalized / weight_plus;

            % compute effective number of particles   
            sum_weight = 0;
            for i = 1 : obj.n
                sum_weight = sum_weight + obj.particle_weight(i) * obj.particle_weight(i);
            end
            Neff = 1 / sum_weight;
            
            if Neff < obj.n / 3
                resample(obj); 
            end
            
            %%%%%%%
            %Calculate the mean and variance from the particles
            meanAndVariance(obj);

        end 
         
        function resample(obj)
            newSamples = zeros(size(obj.particles));
            newWeight = zeros(size(obj.particle_weight));
            W = cumsum(obj.particle_weight);
            r = rand/obj.n;
            count = 1;
            for j = 1:obj.n
                u = r+(j-1)/obj.n;
                while u > W(count)
                    count = count+1;
                end
                newSamples(:,j) = obj.particles(:,count);
                newWeight(j) = 1/obj.n;
            end
            obj.particles = newSamples;
            obj.particle_weight = newWeight;
        end
        
        function meanAndVariance(obj)
            obj.mu = mean(obj.particles, 2); 
            % orientation is a bit more tricky.
            sinSum = 0;
            cosSum = 0;
            for s = 1:obj.n
                cosSum = cosSum + cos(obj.particles(3,s));
                sinSum = sinSum + sin(obj.particles(3,s));
            end
            obj.mu(3) = atan2(sinSum, cosSum);     
            % Compute covariance.
            zeroMean = obj.particles - repmat(obj.mu, 1, obj.n);
            for s = 1:obj.n
                zeroMean(3,s) = wrapTo2Pi(zeroMean(3,s));
            end
            
            obj.Sigma = zeroMean * zeroMean' / obj.n;
        end
    end
end