function[IntValue] = Thrust_Coefficient_Func(varargin)
% Thrust_Coefficient_Func
% 
%
% Output
% IntValue= Valore dell'integrale calcolato per J=.1/pi e Beta075=25/180*pi

if nargin==3
    SigmaValue=varargin{1};
    Beta075Value=varargin{2};
    Jmax=varargin{3};

else
    SigmaValue=0.139;
    Beta075Value=45/180*pi;
    Jmax=1;
end
    
Func=@(x,Lambda,sigma,Beta075)...
[(pi^3.*x.*sqrt(x.^2 + Lambda.^2).*sigma.*(6.*(-atan(Lambda./x) + atan((0.75.*tan(Beta075))./x) +...
(4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma -...
4.*x.*sqrt((4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma).^2./(16.*x.^2) +...
(3.*sqrt(1 + Lambda.^2./x.^2).*sigma.*(-atan(Lambda./x) + atan((0.75.*tan(Beta075))./x)))./x))./(8.*x)) -...
(Lambda.*(0.006 + (0.36.*(0.025.*x - 0.5.*Lambda - 0.375.*sqrt(1 + Lambda.^2./x.^2).*sigma +...
x.*atan(Lambda./x) - 1.*x.*atan((0.75.*tan(Beta075))./x) +...
0.5.*x.*sqrt((4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma).^2/(16.*x.^2) +...
(3.*sqrt(1 + Lambda.^2./x.^2).*sigma.*(-atan(Lambda./x) + atan((0.75.*tan(Beta075))./x)))./x)).^2)/...
x.^2))./x).*(cos((-4.*Lambda - 3.*sqrt(1 + Lambda.^2./x.^2).*sigma +...
4.*x.*sqrt((4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma).^2./(16.*x.^2) +...
(3.*sqrt(1 + Lambda.^2./x.^2).*sigma.*(-atan(Lambda./x) + atan((0.75.*tan(Beta075))./x)))./x))./(8.*x)) -...
(Lambda.*sin((-4.*Lambda - 3.*sqrt(1 + Lambda.^2./x.^2).*sigma +...
4.*x.*sqrt((4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma).^2./(16.*x.^2) +...
(3.*sqrt(1 + Lambda.^2./x.^2).*sigma.*(-atan(Lambda./x) +...
atan((0.75.*tan(Beta075))./x)))./x))./(8.*x)))./x))/8.];



IntValue = integral(@(x)Func(x,.1/pi,SigmaValue,25/180*pi),.1,1);



for i=1:2000
J=linspace(0,Jmax,2000);
y(i)=integral(@(x)Func(x,J(i)/pi,SigmaValue,Beta075Value),.1,1);
end


%figure
hold on
plot(J,y,'-r','LineWidth',3)
xlabel('J')
ylabel('C_t')
grid
hold off

end

