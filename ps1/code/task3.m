close all;
clc;
clear;
clf;

mu_r = 10;
mu_b = 0;
std_rn = 0.5;
std_bn = 0.25;

%%
%3.A

range = repmat(mu_r, 10000 , 1) + randn(10000,1) * std_rn;
bearing = repmat(mu_b, 10000 , 1) + randn(10000,1) * std_bn; 
% x y 
subplot(1,2,2);
X = range .* cos(bearing);
Y = range .* sin(bearing);



plot(X,Y,'g.');
title('Cartesian Frame','fontsize',20,'FontName', 'Times New Roman','Interpreter', 'Latex');

% r-theta(range, bearing, '.');
subplot(1,2,1);
plot(range,bearing,'g.');
title('Sensor Frame','fontsize',20,'FontName', 'Times New Roman','Interpreter', 'Latex');



%%
%3.B Jacobian, Covariance
J = [cos(mu_b) -mu_r * sin(mu_b);
       sin(mu_b)  mu_r * cos(mu_b) ];

covrt = [std_rn^2  0 ; 0  std_bn^2];
mu_rt = [mu_r * cos(mu_b) ; mu_r * sin(mu_b)];
covxy = J * covrt * J.';
mu_xy = J * mu_rt;


%% 3.C
% X-Y
figure;
subplot(1,2,1);
plot(X,Y,'g.');
hold on
%3.C Draw the ellipse (Linearized)
draw_ellipse(mu_xy, covxy , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_xy, covxy , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_xy, covxy , 9, 'r--','linewidth',1.5);

%3.C Sample-based actual
sam_avg_xy = [mean(X) ; mean(Y)];
sam_cov_xy = cov(X,Y); 

draw_ellipse(sam_avg_xy, sam_cov_xy, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 9, 'b--','linewidth',1.5);
title('Cartesian Frame','fontsize',20,'FontName', 'Times New Roman','Interpreter', 'Latex');

%--------------------------------------------------------
%r-theta
subplot(1,2,2);
plot(range,bearing,'g.');
title('Sensor Frame','fontsize',20,'FontName', 'Times New Roman','Interpreter', 'Latex');
hold on
%3.C Draw the ellipse (Linearized)
draw_ellipse(mu_rt, covrt , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt , 9, 'r--','linewidth',1.5);

%3.C Sample-based actual
sam_avg_rt = [mean(range) ; mean(bearing)];
sam_cov_rt = cov(range,bearing); 

draw_ellipse(sam_avg_rt, sam_cov_rt, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 9, 'b--','linewidth',1.5);



%% 3.D

chi2cdf(1,2)
chi2cdf(4,2)
chi2cdf(9,2)

XX = X - 10;
YY = Y - 0;
RR = range - 10;
BB = bearing - 0;

count1_xy = 0;
count2_xy = 0;
count3_xy = 0;
count1_rt = 0;
count2_rt = 0;
count3_rt = 0;

for i = 1 : 10000
Mdis(i) = sqrt([XX(i) ; YY(i)].' * inv(covxy) * [XX(i) ; YY(i)] );
    if Mdis(i) < 1
        count1_xy = count1_xy + 1;
    end
    if Mdis(i) < 2
        count2_xy = count2_xy + 1;
    end
    if Mdis(i) < 3
        count3_xy = count3_xy + 1;
    end
end

for i = 1 : 10000
Mdis(i) = sqrt([RR(i) ; BB(i)].' * inv(covrt) * [RR(i) ; BB(i)] );
    if Mdis(i) < 1
        count1_rt = count1_rt + 1;
    end
    if Mdis(i) < 2
        count2_rt = count2_rt + 1;
    end
    if Mdis(i) < 3
        count3_rt = count3_rt + 1;
    end
end
%% 3.E 
mu_r = 10;
mu_b = 0;
std_rn_tune = 500000; % Should 
std_bn_tune = 0.01; % SHould 





% range bearing
range_tune = repmat(mu_r, 10000 , 1) + randn(10000,1) * std_rn_tune;
bearing_tune = repmat(mu_b, 10000 , 1) + randn(10000,1) * std_bn_tune;
% x y 
X_tune = range .* cos(bearing_tune);
Y_tune = range .* sin(bearing_tune);

XX_tune = X_tune - 10;
YY_tune = Y_tune - 0;
RR_tune = range_tune - 10;
BB_tune = bearing_tune - 0;


count_xy_tune = zeros(3,1);

count_rt_tune = zeros(3,1);

for i = 1 : 10000
Mdis(i) = sqrt([XX(i) ; YY(i)].' * inv(covxy) * [XX(i) ; YY(i)] );
    if Mdis(i) < 1
        count_xy_tune(1) = count_xy_tune(1) + 1;
    end
    if Mdis(i) < 2
        count_xy_tune(2) = count_xy_tune(2) + 1;
    end
    if Mdis(i) < 3
        count_xy_tune(3) = count_xy_tune(3) + 1;
    end
end

for i = 1 : 10000
Mdis(i) = sqrt([RR(i) ; BB(i)].' * inv(covrt) * [RR(i) ; BB(i)] );
    if Mdis(i) < 1
        count_rt_tune(1) = count_rt_tune(1) + 1;
    end
    if Mdis(i) < 2
        count_rt_tune(2) = count_rt_tune(2) + 1;
    end
    if Mdis(i) < 3
        count_rt_tune(3) = count_rt_tune(3) + 1;
    end
end





%% 3.F Using different off-diagonal
% https://en.wikipedia.org/wiki/Covariance_and_correlation


newcov1 = 0.1 * std_rn * std_bn;
newcov2 = 0.5 * std_rn * std_bn;
newcov3 = 0.9 * std_rn * std_bn;

covrt1 = [std_rn^2  newcov1 ; newcov1  std_bn^2];
covxy1 = J * covrt1 * J.';
covrt2 = [std_rn^2  newcov2 ; newcov2  std_bn^2];
covxy2 = J * covrt2 * J.';
covrt3 = [std_rn^2  newcov3 ; newcov3  std_bn^2];
covxy3 = J * covrt3 * J.';


L1 = chol(covrt1,'lower');
L2 = chol(covrt2,'lower');
L3 = chol(covrt3,'lower');
% T1 = sqrt( * inv([0.5 0 ; 0 0.25]) * ) * L1;
% T2 = sqrt( * inv([0.5 0 ; 0 0.25]) * ) * L2;
% T3 = sqrt( * inv([0.5 0 ; 0 0.25]) * ) * L2;


%Z = LX + mu!!!
for i = 1 : 10000 
    d = randn(2,1);
    AA = L1 * d + [10;0];
    range1(i,1) = AA(1,1);
    bearing1(i,1) = AA(2,1);
    %original jacobian transformation
    X1(i,1) = range1(i,1) * cos(bearing1(i,1));
    Y1(i,1) = range1(i,1) * sin(bearing1(i,1));
end

for i = 1 : 10000 
    d = randn(2,1);
    AA = L2 * d + [10;0];
    range2(i,1) = AA(1,1);
    bearing2(i,1) = AA(2,1);

    X2(i,1) = range2(i,1) * cos(bearing2(i,1));
    Y2(i,1) = range2(i,1) * sin(bearing2(i,1));
end
for i = 1 : 10000 
    d = randn(2,1);
    AA = L3 * d + [10;0];
    range3(i,1) = AA(1,1);
    bearing3(i,1) = AA(2,1);
    

    X3(i,1) = range3(i,1) * cos(bearing3(i,1));
    Y3(i,1) = range3(i,1) * sin(bearing3(i,1));
end


figure;
% r-theta(range, bearing, '.');
subplot(2,3,1);
plot(range1,bearing1,'g.');
xlim([4 14])
title('Sensor Frame \rho _{r\theta} = 0.1','fontsize',18,'FontName', 'Times New Roman');
% x y
subplot(2,3,4);
plot(X1,Y1,'g.');
xlim([4 14])
title('Cartesian Frame \rho _{r\theta} = 0.1','fontsize',18,'FontName', 'Times New Roman');

subplot(2,3,2);
plot(range2,bearing2,'g.');
xlim([4 14])
title('Sensor Frame \rho _{r\theta} = 0.5','fontsize',18,'FontName', 'Times New Roman');
% x y
subplot(2,3,5);
plot(X2,Y2,'g.');
xlim([4 14])
title('Cartesian Frame \rho _{r\theta} = 0.5','fontsize',18,'FontName', 'Times New Roman');

subplot(2,3,3);
plot(range3,bearing3,'g.');
xlim([4 14])
title('Sensor Frame \rho _{r\theta} = 0.9','fontsize',18,'FontName', 'Times New Roman');
% x y
subplot(2,3,6);
plot(X3,Y3,'g.');
xlim([4 14])
title('Cartesian Frame \rho _{r\theta} = 0.9','fontsize',18,'FontName', 'Times New Roman');





%3.F(C)--------



%ro = 0.1 --------------------------------------------------------
figure;
subplot(2,3,1);
plot(X1,Y1,'g.');
xlim([5 12])
hold on
%3 Draw the ellipse (Linearized)
draw_ellipse(mu_rt, covxy1 , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covxy1 , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covxy1 , 9, 'r--','linewidth',1.5);
title('Cartesian Frame \rho _{r\theta} = 0.1','fontsize',18,'FontName', 'Times New Roman');
%3 Sample-based actual
sam_avg_xy = [mean(X1) ; mean(Y1)];
sam_cov_xy = cov(X1,Y1); 
draw_ellipse(sam_avg_xy, sam_cov_xy, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 9, 'b--','linewidth',1.5);


subplot(2,3,4);
plot(range1,bearing1,'g.');
xlim([5 12])
title('Sensor Frame \rho _{r\theta} = 0.1','fontsize',18,'FontName', 'Times New Roman');
hold on
%3 Draw the ellipse (r-theta)
draw_ellipse(mu_rt, covrt1 , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt1 , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt1 , 9, 'r--','linewidth',1.5);
%3 Sample-based actual
sam_avg_rt = [mean(range1) ; mean(bearing1)];
sam_cov_rt = cov(range1,bearing1); 
draw_ellipse(sam_avg_rt, sam_cov_rt, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 9, 'b--','linewidth',1.5);
%--------------------------------------------------------


%ro = 0.5 --------------------------------------------------------

subplot(2,3,2);
plot(X2,Y2,'g.');
xlim([5 12])
hold on
%3 Draw the ellipse (linearized)
draw_ellipse(mu_rt, covxy2 , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covxy2 , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covxy2 , 9, 'r--','linewidth',1.5);
title('Cartesian Frame \rho _{r\theta} = 0.5','fontsize',18,'FontName', 'Times New Roman');
%3 Sample-based actual
sam_avg_xy = [mean(X2) ; mean(Y2)];
sam_cov_xy = cov(X2,Y2); 
draw_ellipse(sam_avg_xy, sam_cov_xy, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 9, 'b--','linewidth',1.5);


subplot(2,3,5);
plot(range2,bearing2,'g.');
xlim([5 12])
title('Sensor Frame \rho _{r\theta} = 0.5','fontsize',18,'FontName', 'Times New Roman');
hold on
%3 Draw the ellipse (r-theta)
draw_ellipse(mu_rt, covrt2 , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt2 , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt2 , 9, 'r--','linewidth',1.5);
%3 Sample-based actual
sam_avg_rt = [mean(range2) ; mean(bearing2)];
sam_cov_rt = cov(range2,bearing2); 
draw_ellipse(sam_avg_rt, sam_cov_rt, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 9, 'b--','linewidth',1.5);
%--------------------------------------------------------

%ro = 0.9 --------------------------------------------------------

subplot(2,3,3);
plot(X3,Y3,'g.');
xlim([5 12])
hold on
%3 Draw the ellipse (Linearized)
draw_ellipse(mu_rt, covxy3 , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covxy3 , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covxy3 , 9, 'r--','linewidth',1.5);
title('Cartesian Frame \rho _{r\theta} = 0.9','fontsize',18,'FontName', 'Times New Roman');
%3 Sample-based actual
sam_avg_xy = [mean(X3) ; mean(Y3)];
sam_cov_xy = cov(X3,Y3); 
draw_ellipse(sam_avg_xy, sam_cov_xy, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_xy, sam_cov_xy, 9, 'b--','linewidth',1.5);


subplot(2,3,6);
plot(range3,bearing3,'g.');
xlim([5 12])
title('Sensor Frame \rho _{r\theta} = 0.9','fontsize',18,'FontName', 'Times New Roman');
hold on
%3 Draw the ellipse (Linearized)
draw_ellipse(mu_rt, covrt3 , 1, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt3 , 4, 'r--','linewidth',1.5);
draw_ellipse(mu_rt, covrt3 , 9, 'r--','linewidth',1.5);
%3 Sample-based actual
sam_avg_rt = [mean(range3) ; mean(bearing3)];
sam_cov_rt = cov(range3,bearing3); 
draw_ellipse(sam_avg_rt, sam_cov_rt, 1, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 4, 'b--','linewidth',1.5);
draw_ellipse(sam_avg_rt, sam_cov_rt, 9, 'b--','linewidth',1.5);
%--------------------------------------------------------
%%











function varargout = draw_ellipse(mu, Sigma, k2, varargin)
% h = draw_ellipse(mu, Sigma, K2)
% mu     is the [2x1] [x;y] mean vector
% Sigma  is the [2x2] covariance matrix
% K2     is the Chi-squared 2 DOF variable
% h      is the output plot handle
%
% Draws an ellipse centered at mu with covariance Sigma and
% confidence region k2, i.e.,
% K2 = 1; % 1-sigma
% K2 = 4; % 2-sigma
% K2 = 9; % 3-sigma
% K2 = chi2inv(.50, 2); % 50% probability contour
%
% h = draw_ellipse(mu, Sigma, K2, 'Npoints', N)
% Use the 'Npoints' option to specify the number of points used to define the
% contour, N=20 is the default.
%
% h = draw_ellipse(mu,Sigma,K2,varargin)
% Optional line plot commands can be passed via the varargin syntax, e.g.
% % 1-sigma contour with red dash-dot ellipse and line thickness of 1.5
% draw_ellipse([0;0], eye(2), 1, 'r-.-','linewidth',1.5); 
%
%-------------------------------------------------------
% 2004-11-15    rme    Rewrote to use calculateEllipseXY.m
% 2004-11-21    rme    Added varargin option.
% 2006-05-25    rme    Added optional 'points' argument.
% 2006-06-03    rme    Changed 'points' to 'Npoints' and fixed a bug 
%                      in setting it's default value.
% 2006-06-18    rme    Made plot handle, h, an optional output arg.
% 2009-09-09    rme    Made help description more verbose.

Npoints = 20;
if (nargin > 3)
  for ii=1:length(varargin) % check for 'Npoints', N pair option
    if ischar(varargin{ii}) && strcmpi(varargin{ii},'npoints');
      Npoints = varargin{ii+1};
      varargin = {varargin{1:ii-1},varargin{ii+2:end}};
      break;
    end
  end
end;
    
[x,y] = calculateEllipseXY(mu, Sigma, k2, Npoints);

if nargin > 3
  hp = plot(x,y,varargin{:},'linewidth',1.25); %red
else
  hp = plot(x,y,'linewidth',1.25);  %red
end

if (nargout > 0)
  varargout{1} = hp;
end


end
