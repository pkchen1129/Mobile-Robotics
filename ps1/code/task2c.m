clc; 
clear;
%Task2C


%Two step for now
%% Bel setting
P_land_sensed = 0.8;
P_land_other = 0.4;

bel = zeros(3,10);
bel(1,:) = P_land_other;
bel(1,1) = P_land_sensed;
bel(1,4) = P_land_sensed;
bel(1,7) = P_land_sensed; 
%normalized
bel(1,:) = bel(1,:) / sum(bel(1,:));




%% moves 3 grid cells counterclockwise and detects a landmark
%post_bel = (likelihood * prior)/n
post_bel = zeros(1,10);

for i = 1:10
    if i-3 <= 0    
        shift = i - 3 + 10;
    else
        shift = i - 3; 
    end
    post_bel(i) = bel(1,shift); 
end



for i = 1:10
    if i ==1 || i == 4 || i ==7
        post_bel(i) = post_bel(i) * 0.8;
    else
        post_bel(i) = post_bel(i) * 0.4;
    end
end

post_bel = post_bel/sum(post_bel);
bel(2,:) = post_bel;

%% moves 4 grid cells / Detect none
post_bel = zeros(1,10);
%MAtch the prior
for i = 1:10
    if i-4 <= 0 % Four setps and detect
        shift = i - 4 + 10;
    else
        shift = i - 4;
    end
    post_bel(i) =  bel(2,shift); 
end
for i = 1:10
    if i ==1 || i == 4 || i ==7
        post_bel(i) = post_bel(i) * 0.2; %P(z=1 | x = 1)
    else
        post_bel(i) = post_bel(i) * 0.6; %P(z=1 | x = 0)
    end
end
%Normalized the post bel
post_bel = post_bel / sum(post_bel);
bel(3,:) = post_bel;

