





function [IntValue] = Power_Coefficient_Fun(varargin)
% Power_Coefficient_Fun
%
%
%
% Output
% IntValue= Valore dell'integrale calcolato per J=.2/pi e Beta075=45/180*pi




if nargin==3
    SigmaValue=varargin{1};
    Beta075Value=varargin{2};
    Jmax=varargin{3};
else
    SigmaValue=0.139;
    Beta075Value=45/180*pi;
    Jmax=4;
end


 Func=@(x,Lambda,Sigma,Beta075)...
[(pi^4.*x.^2.*sqrt(x.^2 +Lambda.^2).*Sigma.*(0.006 + (3.*Lambda.* ...
(4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*Sigma - 8.*x.*atan(Lambda./x) + ...
8.*x.*atan((0.75.*tan(Beta075))./x) - 4.*x.*sqrt((4.*Lambda + ...
3.*sqrt(1 +Lambda.^2./x.^2).*Sigma).^2./(16.*x.^2) + ...
(3.*sqrt(1 +Lambda.^2./x.^2).*Sigma.* (-atan(Lambda./x) +...
atan((0.75.*tan(Beta075))./x)))./x)))./(4.*x.^2) + ...
(0.36.*(0.025.*x - 0.5.*Lambda - ...
0.375.*sqrt(1 +Lambda.^2./x.^2).*Sigma + x.*atan(Lambda./x) - ...
1.*x.*atan((0.75.*tan(Beta075))./x) + 0.5.*x.*sqrt((4.*Lambda + ...
3.*sqrt(1 + Lambda.^2./x.^2).*Sigma).^2./(16.*x.^2) + ...
(3.*sqrt(1 +Lambda.^2./x.^2).*Sigma.*(-atan(Lambda./x) + ...
atan((0.75.*tan(Beta075))./x)))./x)).^2)./x.^2).* ...
(cos((-4.*Lambda - 3.*sqrt(1 + Lambda.^2./x.^2).*Sigma + 4.*x.*sqrt((4.*Lambda +... 
3.*sqrt(1 + Lambda.^2./x.^2).*Sigma).^2./(16.*x.^2) + ...
(3*sqrt(1 + Lambda.^2./x.^2).*Sigma.*(-atan(Lambda./x) + ...
atan((0.75.*tan(Beta075))./x)))./x))./(8.*x)) - ...
(Lambda.*sin((-4.*Lambda - 3.*sqrt(1 + Lambda.^2./x.^2).*Sigma + 4.*x.* ...
sqrt((4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*Sigma).^2./(16.*x.^2) + ...
(3.*sqrt(1 + Lambda.^2./x.^2).*Sigma.*(-atan(Lambda./x) + ...
atan((0.75.*tan(Beta075))./x)))./x))./(8.*x)))./x)./8)];



IntValue = integral(@(x)Func(x,.2/pi,SigmaValue,45/180*pi),.1,1);



for i=1:2000
J=linspace(0,Jmax,2000);
y(i)=integral(@(x)Func(x,J(i)/pi,SigmaValue,Beta075Value),.1,1);
end


%figure
hold on
plot(J,y,'-b','LineWidth',3)
xlabel('J')
ylabel('C_p')
grid
hold off

end