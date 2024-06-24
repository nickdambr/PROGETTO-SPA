% Tabella dell'aria tipo
H=h;
z = [                                                                      % z = quota [m]
    0
    500
    1000
    1500
    2000
    2500
    3000
    3500
    4000
    4500
    5000
    5500
    6000
    6500
    7000
    7500
    8000
    8500
    9000
    9500
    10000
    10500
    11000
    11500
    12000
    12500
    13000
    13500
    13500
    14000
    14500
    15000
    15500
    16000
    16500
    17000
    17500
    18000
    18500
    19000
    19500
    ];

P = [                                                                      % P = Pressione [Pa]
    101325.25
    95450.10
    89853.89
    84526.38
    79457.57
    74637.76
    70057.47
    65707.49
    61578.85
    57662.82
    53950.94
    50434.95
    47106.86
    43958.87
    40983.43
    38173.22
    35521.11
    33020.22
    30663.86
    28445.55
    26359.01
    24398.17
    22557.15
    20847.42
    19266.78
    17805.99
    16455.95
    15208.53
    14055.19
    12989.53
    12004.67
    11094.48
    10253.31
    9475.91
    8757.45
    8093.46
    7479.82
    6912.71
    6388.59
    5904.21
    5456.55
    ];

T = [                                                                      % T = Temperatura [K]
    288.16
    284.91
    281.66
    278.41
    275.16
    271.91
    268.66
    265.41
    262.16
    258.91
    255.66
    252.41
    249.16
    245.91
    242.66
    239.41
    236.16
    232.91
    229.66
    226.41
    223.16
    219.91
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    216.66
    ];

rho0 = 1.225;                                                              % [kg/m^3]
beta = 10^-4;                                                              % [m^-1]
rhoa = rho0*(exp(-beta*h));

if ismember(h,z)
   
   ind = find(z == h);
   Pa  = P(ind);
   Ta  = T(ind);

else
    
    l = length(z);
   
    for i = 1 : l
       
       if z(i) > h
           imax = i;
           imin = i-1;
           break
       
       end
       
     end
   
   Ta = T(imin) + ( ((h-z(imin))*(T(imax)-T(imin))) / (z(imax)-z(imin)) ); % Interpolazione tra z e T
   Pa = P(imin) + ( ((h-z(imin))*(P(imax)-P(imin))) / (z(imax)-z(imin)) ); % Interpolazione tra z e P

end
Pa=Pa*10^(-3)
%fprintf('Ta: %f \n',Ta)
%fprintf('Pa: %f \n', Pa)
%fprintf('rhoa: %f \n', rhoa)
