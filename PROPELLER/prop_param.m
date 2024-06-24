%Propeller parameters

%input
v_suono=a;
v = V0;
ngiri = n; %al secondo
Diam = D;
corda = c;
Npale = B;


omega=(2*pi*ngiri);
R=Diam/2;

%output
Mtip=(omega*R)/v_suono; %mach al tip
J=v/(ngiri*Diam); %rapporto di avanzamento
sigma=(Npale*corda)/(pi*R); %rapporto di solidità