classdef PF < handle
    properties
        Q;                  %Sensor Noise
        n;                  %Number of Particles
        particles;          %Pose of particle
        particle_weight;    %Particle Weight
        mu;
        Sigma;
    end
    
    methods
        function obj = PF(init)
            % motion model
            % motion noise covariance
            % measurement noise covariance
            obj.Q;
            % PF parameters
            obj.mu = init.mu;
            obj.Sigma = init.Sigma;
            obj.n = init.n; %Number of particle
            obj.particles = init.particles; %(2,n)
            obj.particle_weight = init.particle_weight;
        end
        function prediction(obj)
            for j = 1:obj.n
                obj.particles(:,j) = [obj.particles(1,j) + randn(1);
                                      obj.particles(2,j) + randn(1)];
            end
            
        end
        
        
        function correction(obj, z)
            %data type of z is like z = [r1, r2, r3];
            landmark_x = [-100; 60; -20];
            landmark_y = [  30; 20; -40];
            
%             u_noise_std = chol(obj.M(u), 'lower');
            weight = zeros(obj.n,1);
            for j = 1:obj.n
                %Propagate Particles through Motion Model
%                 obj.particles(:,j) = obj.gfun(obj.particles(:,j),u_noise_std*randn(3,1)+u);
                %Calculate measurement and difference in measurements
%                 z_hat = obj.hfun(landmark_x, landmark_y, obj.particles(:,j));
                z_hat = [ sqrt((landmark_x(1) - obj.particles(1,j))^2 + (landmark_y(1) - obj.particles(2,j))^2) ;...
                          sqrt((landmark_x(2) - obj.particles(1,j))^2 + (landmark_y(2) - obj.particles(2,j))^2) ;...
                          sqrt((landmark_x(3) - obj.particles(1,j))^2 + (landmark_y(3) - obj.particles(2,j))^2) ];
                
                diff = [...
                        z(1) - z_hat(1);
                        z(2) - z_hat(2);
                        z(3) - z_hat(3) ];
                %Use mvnpdf to get weight probability of difference in measurement
                if z(3) < 25
                    obj.Q = diag([1 1 1]);
                else
                    obj.Q = diag([z(1).^(1/3) z(2).^(1/3) z(3).^(1/3)]);
                end
                
                weight(j) = mvnpdf(diff, 0, 2 * obj.Q);
            end
            %Update Weights
            
            obj.particle_weight = obj.particle_weight.*weight;
            obj.particle_weight = obj.particle_weight./sum(obj.particle_weight);
            Neff = 1 / sum(obj.particle_weight.^2);
            if Neff < obj.n /5
                obj.resample();
            end
            obj.meanAndVariance();
        end 
         
        
        %% Function resample and meanandvariance
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
            % Compute covariance.
            zeroMean = obj.particles - repmat(obj.mu, 1, obj.n);
            
            obj.Sigma = zeroMean * zeroMean' / obj.n;
        end
    end
end