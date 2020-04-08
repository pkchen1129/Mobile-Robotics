function filter = filter_initialization(initialStateMean,initialStateCov, filter_name)
switch filter_name
    case "PF"
        init.mu =  initialStateMean;
        init.Sigma = 0.01 * eye(2);
        init.n = 200;
        init.particles = zeros(2, init.n);
        init.particle_weight = zeros(init.n, 1);
        L = chol(init.Sigma,'lower');
        for i = 1:init.n
            init.particles(:,i) = L*randn(2,1) + init.mu(:,1);
            init.particle_weight(i) = 1/init.n;
        end
        filter = PF(init);
        

end
end