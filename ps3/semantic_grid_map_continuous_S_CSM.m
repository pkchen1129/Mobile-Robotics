classdef semantic_grid_map_continuous_S_CSM < handle
    properties
        % map dimensions
        range_x = [-15, 20];
        range_y = [-25, 10];
        % sensor parameters
        z_max = 30;                 % max range in meters
        n_beams = 133;              % number of beams
        % grid map paremeters
        grid_size = 0.5;           % 0.135 , 0.27 , 0.5
        alpha = 2 * 0.5;           % 2 * grid_size
        beta = 2 * pi/133;          % 2 * pi/n_beams
        nn = 16;                    % number of nearest neighbor search
        map;                        % map!
        pose;                       % pose data
        scan;                       % laser scan data
        m_i = [];                   % cell i
        % semantic
        num_classes = 6;
        
        l = 0.2;
        sigma = 0.1;
        % -----------------------------------------------
        % To Do: 
        % prior initialization
        alpha0 = 0.001 ;
        
        % -----------------------------------------------
    end
    
    methods
        function obj = semantic_grid_map_continuous_S_CSM(pose, scan)
            % class constructor
            % construct map points, i.e., grid centroids.
            x = obj.range_x(1):obj.grid_size:obj.range_x(2);
            y = obj.range_y(1):obj.grid_size:obj.range_y(2);
            [X,Y] = meshgrid(x,y);
            t = [X(:), Y(:)];
            % a simple KDtree data structure for map coordinates.
            obj.map.occMap = KDTreeSearcher(t);
            obj.map.size = size(t,1);
            
            % -----------------------------------------------
            % To Do: 
            % map parameter initialization such as map.alpha
            
            obj.map.alpha = obj.alpha0 * ones(obj.map.size , obj.num_classes + 1);
            obj.map.mean = zeros(obj.map.size , obj.num_classes + 1);
            obj.map.variance = zeros(obj.map.size , 1);
            % -----------------------------------------------
            
            % set robot pose and laser scan data
            obj.pose = pose;
            obj.pose.mdl = KDTreeSearcher([pose.x, pose.y]);
            obj.scan = scan;
        end
        
        function build_ogm(obj)
            % build occupancy grid map using the binary Bayes filter.
            % we first loop over all map cells, then for each cell, we find
            % N nearest neighbor poses to build the map. Note that this is
            % more efficient than looping over all poses and all map cells
            % for each pose which should be the case in online
            % (incremental) data processing.
            for i = 1:obj.map.size
                m = obj.map.occMap.X(i,:);
                idxs = knnsearch(obj.pose.mdl, m, 'K', obj.nn);
                if ~isempty(idxs)
                    for k = idxs
                        % pose k
                        pose_k = [obj.pose.x(k),obj.pose.y(k), obj.pose.h(k)];
                        if obj.is_in_perceptual_field(m, pose_k)
                            % laser scan at kth state; convert from
                            % cartesian to polar coordinates
                            [bearing, range] = cart2pol(obj.scan{k}(1,:), obj.scan{k}(2,:));
                            z = [range' bearing' obj.scan{k}(3,:)'];
                            
                            % -----------------------------------------------
                            % To Do: 
                            % update the semantic sensor model in cell i
                            obj.continuous_S_CSM(z, i);
                            
                            % -----------------------------------------------
                        end
                    end
                end
                
                % -----------------------------------------------
                % To Do: 
                % update mean and variance for each cell i
                N = sum(obj.map.alpha(i,:),2);
                obj.map.mean(i,:) = obj.map.alpha(i,:) / N; %1*7
                max_mean = max(obj.map.mean(i,:));  %1*1
                obj.map.variance(i) = max_mean * (1 - max_mean) / (N + 1);
                % -----------------------------------------------
                
            end
        end
        
        function inside = is_in_perceptual_field(obj, m, p)
            % check if the map cell m is within the perception field of the
            % robot located at pose p.
            inside = false;
            d = m - p(1:2);
            obj.m_i.range = sqrt(sum(d.^2));
            obj.m_i.phi = wrapToPi(atan2(d(2),d(1)) - p(3));
            % check if the range is within the feasible interval
            if (0 < obj.m_i.range) && (obj.m_i.range < obj.z_max)
                % here sensor covers -pi to pi!
                if (-pi < obj.m_i.phi) && (obj.m_i.phi < pi)
                    inside = true;
                end
            end
        end
        
        function continuous_S_CSM(obj, z, i)
            % -----------------------------------------------
            % To Do: 
            % implement the continuous semantic counting sensor model
            bearing_diff = abs(wrapToPi(z(:,2) - obj.m_i.phi));
            [~, k] = min(bearing_diff);
            
            x_star = [obj.m_i.range * cos(obj.m_i.phi) , obj.m_i.range * sin(obj.m_i.phi)]; %1*2
            x_i = [ z(k,1) .* cos(z(k,2)) , z(k,1) .* sin(z(k,2)) ]; %1*2
            
            %Set d_alpha
            d_alpha = sqrt(sum((x_i - x_star).^2));
            
            %Set condition for projection
            PQlen = sqrt(sum(x_i.^2)); %1*1
            projection_len = (x_star * x_i' / PQlen ); %length of vector P-x_free
            projection = projection_len * (x_i/norm(x_i)); %vector P-x_free(1,2)
            
            if projection_len <= PQlen & projection_len >= 0
                x_free = projection;
                 %1*1
                d_beta = sqrt(sum((x_free - x_star).^2));
                
                if (d_alpha < obj.l)
                    kernel = obj.sigma *...
                    ( (2+cos(2*pi*d_alpha/obj.l))/3 * (1 - d_alpha/obj.l) + 1/(2*pi)*sin(2*pi*d_alpha/obj.l));
                    obj.map.alpha(i , z(k,3)) = obj.map.alpha(i , z(k,3)) + kernel;
                   
                    
                elseif (d_beta < obj.l)
                    kernel = obj.sigma *...
                    ( (2+cos(2*pi*d_beta/obj.l))/3 * (1 - d_beta/obj.l) + 1/(2*pi)*sin(2*pi*d_beta/obj.l));
                    obj.map.alpha(i , 7) = obj.map.alpha(i , 7) + kernel;
                    
                end

            elseif projection_len > PQlen
                if (d_alpha < obj.l)
                    kernel = obj.sigma *...
                    ( (2+cos(2*pi*d_alpha/obj.l))/3 * (1 - d_alpha/obj.l) + 1/(2*pi)*sin(2*pi*d_alpha/obj.l));
                    obj.map.alpha(i , z(k,3)) = obj.map.alpha(i , z(k,3)) + kernel;
                
                end
            elseif projection_len < 0
                %do nothing
            end
            
            
            
            
            
            
            
            
            
            % -----------------------------------------------
        end
    end
end