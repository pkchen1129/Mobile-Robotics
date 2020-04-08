function varargout = run_randomwalk(filter_name)
% Using filter_initialization.m to initialize the value and Pf.m to spcify
% the prediction and correction
%
addpath([cd, filesep, 'lib'])

cla;
%--------------------------------------------------------------
% Initializations
%--------------------------------------------------------------
initialStateMean = [100 -10]';
initialStateCov = 0.01 * eye(2);
persistent numSteps;
numSteps = 200;
z = load('measurements.mat').Z;

results = zeros(2,numSteps);
% sys = system_initialization();
filter = filter_initialization(initialStateMean, initialStateCov, filter_name);
z = load('measurements.mat').Z;

for t = 1:numSteps
    %=================================================
    % data available to your filter at this time step
    %=================================================
    
    observation = z(t,1:3)';   % [x1, x2, x3]' noisy observation
    
    switch filter_name
        case "PF"
%             hp = plot(filter.particles(1,:), filter.particles(2,:),'.','Color', [[0.2980 .6 0], .25]);
            filter.prediction();
            filter.correction(observation);
%             results(:,t) = mahalanobis(filter,data(t,10:12));
            results(:,t)= filter.mu;
    end

end

X_start = [-100,60, -20];
Y_start = [ 30 ,20, -40];

plot(X_start,Y_start, '*');
hold on
plot(results(1,:),results(2,:));
xlim([-120 120])
ylim([-120 120])
grid on
axis equal
varargout{1} = results;

end

