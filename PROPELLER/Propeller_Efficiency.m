function [outputArg1,outputArg2] = Propeller_Efficiency(varargin)
% Propeller_Efficiency


if nargin==3
    SigmaValue=varargin{1};
    Beta075Value=varargin{2};
    Jmax=varargin{3};
    Color='g';
    
elseif nargin==4
    SigmaValue=varargin{1};
    Beta075Value=varargin{2};
    Jmax=varargin{3};
    Color=varargin{4};
    
else
    SigmaValue=0.139;
    Beta075Value=30/180*pi;
    Jmax=4;
    Color='g';
end

C_t_Func=@(x,Lambda,sigma,Beta075)...
[(pi^3.*x.*sqrt(x.^2 + Lambda.^2).*sigma.*(6.*(-atan(Lambda./x) + atan((0.75.*tan(Beta075))./x) +...
(4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma -...
4.*x.*sqrt((4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*sigma).^2./(16.*x.^2) +...
(3.*sqrt(1 + Lambda.^2./x.^2).*sigma.*(-atan(Lambda./x) + atan((0.75.*tan(Beta075))./x)))./x))./(8.*x)) -...
(Lambda.*(0.006 + (0.36.*(0.024999999999999998.*x - 0.5.*Lambda - 0.375.*sqrt(1 + Lambda.^2./x.^2).*sigma +...
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

C_p_Func=@(x,Lambda,Sigma,Beta075)...
[(pi^4.*x.^2.*sqrt(x.^2 +Lambda.^2).*Sigma.*(0.006 + (3.*Lambda.* ...
(4.*Lambda + 3.*sqrt(1 + Lambda.^2./x.^2).*Sigma - 8.*x.*atan(Lambda./x) + ...
8.*x.*atan((0.75.*tan(Beta075))./x) - 4.*x.*sqrt((4.*Lambda + ...
3.*sqrt(1 +Lambda.^2./x.^2).*Sigma).^2./(16.*x.^2) + ...
(3.*sqrt(1 +Lambda.^2./x.^2).*Sigma.* (-atan(Lambda./x) +...
atan((0.75.*tan(Beta075))./x)))./x)))./(4.*x.^2) + ...
(0.36.*(0.024999999999999998.*x - 0.5.*Lambda - ...
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



eta_p_Func=zeros(1,2000);
for i=1:2000
    J=linspace(0,Jmax,2000);
    eta_p_Func(i)=((integral(@(x)C_t_Func(x,J(i)/pi,SigmaValue,Beta075Value),.1,1)).*J(i))./...
        (integral(@(x)C_p_Func(x,J(i)/pi,SigmaValue,Beta075Value),.1,1));
    if i>1
        if eta_p_Func(i)<=0
            eta_p_Func(i)=0;
            break
        end
    end
end

if Color=='g'
hold on
plot(J,eta_p_Func,'-g','LineWidth',3)
xlabel('J')
ylabel('\eta_p')
grid
hold off
elseif Color=='r'
    hold on
plot(J,eta_p_Func,'-r','LineWidth',3)
xlabel('J')
ylabel('\eta_p')
grid
hold off
elseif Color=='b'
    hold on
plot(J,eta_p_Func,'-b','LineWidth',3)
xlabel('J')
ylabel('\eta_p')
grid
hold off

end
end

