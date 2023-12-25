clc;clear;close;
G=1;
C_b=100;
l = [ 8 10 12 14 16 18 20 ];

R = 0:0.1:15;
% lambda = 50;
lambda=100;
for i=1:length(l)
    sigma2= v_3D(l(i),R,lambda,C_b,G);
    plot(R, sigma2)
    ylim([-5 5])
    %axis([0 20 -5 5])    
    hold on
end
    plot(R, zeros(size(R)),'b--')
title('3D in vivo \lambda = 100')
legend('l=8','l=10','l=12','l=14','l=16','l=18','l=20','y=0')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off;

R = 0:0.01:1;
% lambda = 50;
lambda=100;
for i=1:length(l)
    sigma2= v_3D(l(i),R,lambda,C_b,G);
    plot(R, sigma2)
% axis([0 20 -5 5])    
    hold on
end
title('3D in vivo \lambda = 100')
legend('l=8','l=10','l=12','l=14','l=16','l=18','l=20')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off;

R = 0:0.1:50;
% lambda = 1;
lambda=1;
for i=1:length(l)
    sigma2= v_3D(l(i),R,lambda,C_b,G);
    plot(R, sigma2)
    %axis([xmin xmax ymin ymax])
    hold on
end
title('3D in vivo \lambda = 1')
legend('l=8','l=10','l=12','l=14','l=16','l=18','l=20')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off;

R = 0:0.1:50;
% lambda=0.8;
lambda=0.8;
for i=1:length(l)
    sigma2= v_3D(l(i),R,lambda,C_b,G);
    plot(R, sigma2)
    %axis([xmin xmax ymin ymax])
    hold on
end
title('3D in vivo \lambda = 0.8')
legend('l=8','l=10','l=12','l=14','l=16','l=18','l=20')
legend('location','best')
xlabel('R')
ylabel('\delta^{-1}d\delta/dt')
hold off;


%%
function sigma= v_3D(L,R,lambda,C_b,G)
    sigma_A = -(G.*((( -2.*R.*((3.*R+2).*lambda.^(3/2)+sqrt(lambda).*(R-2) ).*(L+1/2).*besselk(L+1/2, R) ...
        -2.*(R.*(R-L./2-1).*lambda.^(3/2)-sqrt(lambda).*(L+2).*(R-2)./2 ).*(R+1).*besselk(L-1/2, R)) ...
        .* (cosh(sqrt(lambda).*R)) .^3 -2.* sinh(sqrt(lambda).*R).* (R.*(L+1/2) ...
        .*((lambda.^2+3.*lambda).*R+lambda.^2-1) .*besselk(L+1/2, R)+ ((lambda.^2+lambda).*R.^2-2.*lambda.*(L+2).*R ...
        +(lambda+1).*(L+2)).*(R+1).*besselk(L-1/2, R)./2).*(cosh(sqrt(lambda).*R)).^2+(2.*R.*((2.*R+2).*lambda.^(3/2) ...
        +sqrt(lambda).*(R-2)).*(L+1/2).*besselk(L+1/2, R)+2.*(R.^2.*lambda.^(3/2)-sqrt(lambda).*(L+2).*(R-2)./2).*(R+1) ...
        .*besselk(L-1/2, R)).*cosh(sqrt(lambda).*R)+2.*sinh(sqrt(lambda).*R).*(R.*(L+1/2).*(lambda.*R+lambda-1) ...
        .*besselk(L+1/2, R)+besselk(L-1/2, R).*(R+1).*(R.^2.*lambda+L+2)./2)).*besseli(L+1/2, sqrt(lambda).*R) ...
        +besseli(L-1/2, sqrt(lambda).*R).*(R.^2+(L+2).*R+L+2).*(lambda.*(-2+(lambda+1).*R) ...
        .*(cosh(sqrt(lambda).*R)) .^3+2.* sinh(sqrt(lambda).*R).*((R-1/2).*lambda.^(3/2)-sqrt(lambda) ./2 ) ...
        .*(cosh(sqrt(lambda).*R)) .^2-lambda.*(R-2).*cosh(sqrt(lambda).*R)+sinh(sqrt(lambda).*R).*sqrt(lambda)) .*besselk(L+1/2, R)).*C_b);
 
     sigma_B = (sqrt(lambda).*R.^3.*( besselk(L+1/2, R).*(cosh(sqrt(lambda).*R).*lambda+sinh(sqrt(lambda).*R).*sqrt(lambda)).*besseli(L-1/2,sqrt(lambda).*R) ...
    + besseli(L+1/2,sqrt(lambda).*R).*besselk(L-1/2, R).*( sqrt(lambda).*cosh(sqrt(lambda).*R) +sinh(sqrt(lambda).*R))) ...
    .*(cosh(sqrt(lambda).*R).*lambda+sinh(sqrt(lambda).*R).*sqrt(lambda)).*( sqrt(lambda).*cosh(sqrt(lambda).*R) +sinh(sqrt(lambda).*R)));
    
    sigma = sigma_A ./ sigma_B;
end