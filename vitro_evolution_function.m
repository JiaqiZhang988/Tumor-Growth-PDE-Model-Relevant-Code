clc;clear;close;

R = 1:0.1:50; % 注意一下这个步长的选择
G0 = 1;
c_B = 10;
%% 
% 

%%fix m, change l
%m1 = 0.5
%m1=13;
l = [3 8 10 12 14 16 18];

% lambda = 5;
lambda = 5;

for i=1:length(l)
    sigma1= vitro_3D(l(i),R,l(i),lambda,c_B,G0);
    plot(R, sigma1)
    hold on
end

title('3D in vitro \lambda = 5')
legend('l=3','l=8','l=10','l=12','l=14','l=16','l=18')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off

% lambda = 1;
lambda = 1;

for i=1:length(l)
    sigma1= vitro_3D(l(i),R,l(i),lambda,c_B,G0);
    plot(R, sigma1)
    hold on
end

title('3D in vitro \lambda = 1')
legend('l=3','l=8','l=10','l=12','l=14','l=16','l=18')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off

% lambda = 0.8;
lambda = 0.8;

for i=1:length(l)
    sigma1= vitro_3D(l(i),R,l(i),lambda,c_B,G0);
    plot(R, sigma1)
    hold on
end

title('3D in vitro \lambda = 0.8')
legend('l=3','l=8','l=10','l=12','l=14','l=16','l=18')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off

% lambda = 0.4;
lambda = 0.4;

for i=1:length(l)
    sigma1= vitro_3D(l(i),R,l(i),lambda,c_B,G0);
    plot(R, sigma1)
    hold on
end

title('3D in vitro \lambda = 0.4')
legend('l=3','l=8','l=10','l=12','l=14','l=16','l=18')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off

% lambda = 0.1;
lambda = 20;

for i=1:length(l)
    sigma1= vitro_3D(l(i),R,l(i),lambda,c_B,G0);
    plot(R, sigma1)
    hold on
end

title('3D in vitro \lambda = 20')
legend('l=3','l=8','l=10','l=12','l=14','l=16','l=18','l=20')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off
%% 
% 


% %%fix l change m
% l=13;
% m1=[8 10 12 14 16 18];
% 
% 
% lambda = 100;
% %c_B = 10;
% 
% for i=1:length(m1)
%     sigma1= vitro_3D(l,R,m1(i),lambda,c_B,G0);
%     plot(R, sigma1)
%     hold on
% end
% 
% title('3D in vitro \lambda = 100')
% legend('m1=8','m1=10','m1=12','m1=14','m1=16','m1=18')
% legend('location','best')
% xlabel('R')
% ylabel('\delta^{-1}d\delta/dt')
% hold off
% 
% 
% lambda = 50;
% %c_B = 10;
% 
% for i=1:length(m1)
%     sigma1= vitro_3D(l,R,m1(i),lambda,c_B,G0);
%     plot(R, sigma1)
%     hold on
% end
% 
% title('3D in vitro \lambda = 50')
% legend('m1=8','m1=10','m1=12','m1=14','m1=16','m1=18')
% legend('location','best')
% xlabel('R')
% ylabel('\delta^{-1}d\delta/dt')
% hold off
% 
% lambda = 2;
% 
% for i=1:length(m1)
%     sigma1= vitro_3D(l,R,m1(i),lambda,c_B,G0);
%     plot(R, sigma1)
%     hold on
% end
% 
% title('3D in vitro \lambda = 2')
% legend('m1=8','m1=10','m1=12','m1=14','m1=16','m1=18')
% legend('location','best')
% xlabel('R')
% ylabel('\delta^{-1}d\delta/dt')
% hold off
% 
% 
% lambda = 0.8;
% %c_B = 10;
% 
% for i=1:length(m1)
%     sigma1= vitro_3D(l,R,m1(i),lambda,c_B,G0);
%     plot(R, sigma1)
%     hold on
% end
% 
% title('3D in vitro \lambda = 0.8')
% legend('m1=8','m1=10','m1=12','m1=14','m1=16','m1=18')
% legend('location','best')
% xlabel('R')
% ylabel('\delta^{-1}d\delta/dt')
% hold off
% 
% lambda = 0.4;
% %c_B = 10;
% 
% for i=1:length(m1)
%     sigma1= vitro_3D(l,R,m1(i),lambda,c_B,G0);
%     plot(R, sigma1)
% 
% 
% 
%     hold on
% end
% 
% title('3D in vitro \lambda = 0.4')
% legend('m1=8','m1=10','m1=12','m1=14','m1=16','m1=18')
% legend('location','best')
% xlabel('R')
% ylabel('\delta^{-1}d\delta/dt')
% hold off

%%
function sigma1= vitro_3D(l,R,m,lambda,c_B,G0)
    A = (R.^(-1-m).* sinh(sqrt(lambda).*R) .* R.^(m+5/2) .* m .* sqrt(lambda))./2 ;
    B = (sqrt(lambda).*sqrt(R) + (lambda.^(3/2).* R.^(5/2))./2 ).* sinh(sqrt(lambda).*R) - cosh(sqrt(lambda).*R) .*lambda .* R.^(3/2);
    C1 =  -lambda .* (l+1) .* besseli(l+1/2,sqrt(lambda).*R);
    C2 =  besseli(l-1/2,sqrt(lambda).*R) .* lambda.^(3/2).*R  ;
    C =  ( C1 + C2).* R.^(5/2) .* cosh(sqrt(lambda).*R) ;
    % C =  ( C1 + C2).* R.^(5/2) .* cosh(sqrt(lambda).*R) ;
    D =  lambda^(3/2) .* R.^(7/2) .* sinh(sqrt(lambda).*R) .* besseli(l+1/2,sqrt(lambda).*R)  ;
    sigma1 =2.*c_B.*G0.* ( (A+B.*R) .* besseli(l+1/2,sqrt(lambda).*R)   - C./2 )./D;
end