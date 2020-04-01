% close all
% cla
% % 2D Raw Data
% xx = inputINTELg2o(1:1227, 3);
% yy = inputINTELg2o(1:1227, 4);
% x = inputINTELg2o(1229:2711, 4);
% y = inputINTELg2o(1229:2711, 5);
% % 2D Batch Optimize
% xxx = out(1:1228,3);
% yyy = out(1:1228,4);
% % 2D Incremental
% xxxx = outisam2(1:1228, 3);
% yyyy = outisam2(1:1228, 4);
% 
% plot(xx,yy,'b.',xxx,yyy,'r.'); %
% legend('Unoptimized','Optimized')
% % plot(xxx,yyy,'color','r');
% 
% axis equal



%-----------------------------------
% 3D Raw Data
xx = parkinggarage(1:1661, 3);
yy = parkinggarage(1:1661, 4);
zz = parkinggarage(1:1661, 5);

% 3D Batch Optimize
xxx = outbatch3D(1:1661, 3);
yyy = outbatch3D(1:1661, 4);
zzz = outbatch3D(1:1661, 5);

% 3D Incremental
xxxx = outisam23D(1:1661, 3);
yyyy = outisam23D(1:1661, 4);
zzzz = outisam23D(1:1661, 5);
% plot3();
axis equal
plot3(xx,yy,zz);
hold on
plot3(xxx,yyy,zzz);


legend('Unoptimized','Optimized')



