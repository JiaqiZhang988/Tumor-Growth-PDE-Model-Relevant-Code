R=0:0.1:20;

plot(R,S50(R,0.5),R,S50(R,0.8),R,S50(R,1),R,S50(R,10),R,S50(R,20))
legend('\lambda=0.5','\lambda=0.8','\lambda=1','\lambda=10','\lambda=20')
legend('location','best')
ylim([0 1.4])
xlabel('R')
ylabel('dR/dt')
title('3D in vivo')

plot(R,S41(R,0.5),R,S41(R,0.8),R,S41(R,1),R,S41(R,10),R,S41(R,20))
legend('\lambda=0.5','\lambda=0.8','\lambda=1','\lambda=10','\lambda=20')
legend('location','best')
%ylim([0 1])
xlabel('R')
ylabel('dR/dt')
title('3D in vitro')
%%
function s= S50(R,lambda)
    s1=(R+1).*(cosh(sqrt(lambda).*R).*R.*lambda-sinh(sqrt(lambda).*R).*sqrt(lambda));
    s2= R.^2.*lambda.*(cosh(sqrt(lambda).*R).*lambda + sinh(sqrt(lambda).*R).*sqrt(lambda)) ;
    s=s1./s2;
end

function s= S41(R,lambda)
    s=cosh(sqrt(lambda).*R)./(sinh(sqrt(lambda).*R).*sqrt(lambda))-1./(R.*lambda);
end