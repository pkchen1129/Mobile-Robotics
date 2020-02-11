clc;
clear;
sigma_x1 = 1;
sigma_x2 = 0.8;
%prior value
mu_o = 0;
sigma_o = sqrt(1000);


s1 = [10.6715 8.7925 10.7172 11.6302 10.4889 11.0347 10.7269 9.6966 10.2939 9.2127];
s2 = [10.7107 9.0823 9.1449 9.3524 10.2602];



for i = 1:10
    mu_o = (sigma_o^2 * s1(i) + sigma_x1^2 * mu_o ) /(sigma_o^2 + sigma_x1^2)
    sigma_o = sigma_x1^2 * sigma_o^2 / (sigma_o^2 + sigma_x1^2)
end




% 
% for i = 1:5
%     mu_o = (sigma_o^2 * s2(i) + sigma_x2^2 * mu_o ) /(sigma_o^2 + sigma_x2^2)
%     sigma_o = sigma_x2^2 * sigma_o^2 / (sigma_o^2 + sigma_x2^2) 
% end
