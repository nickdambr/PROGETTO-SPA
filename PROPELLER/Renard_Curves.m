% Renard_Curves


% Grafico 1

Sigma=0.139;
Beta075=25/180*pi;
Jmax=4;

figure(1)
Thrust_Coefficient_Func(Sigma,Beta075,Jmax);
Power_Coefficient_Fun(Sigma,Beta075,Jmax);
Propeller_Efficiency(Sigma,Beta075,Jmax)

xlabel('J')
ylabel('C_t (Red)-C_p (Blue)_\eta_p (Green)')

axis([0 1.6 0 1])


% Grafico 2

figure(2)

Beta075=40/180*pi;
Propeller_Efficiency(Sigma,Beta075,Jmax,'r')

Beta075=10/180*pi;
Propeller_Efficiency(Sigma,Beta075,Jmax,'b')
