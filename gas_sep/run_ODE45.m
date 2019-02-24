function [P] = run_ODE45(  )
 
zspan=[0 5];

W0=[0.042378;0.159423;0.4158]; %based on concentraion of air at 5 bar, and Q = [L/hr] but # based on 0.01ft/s through 2450-fibre bundle

[z,W_s]=ode45(@(z,W) mem_gas_sep(z,W),zspan,W0);

P = 2450.*W_s(:,1).*W_s(:,3); %o2
P1 = 2450.*W_s(:,2).*W_s(:,3); %n2
P2 = 0.08314*298.*(W_s(:,1)+W_s(:,2)); %delta P
P3 = 2450*(W_s(:,3)); %volumetric flow rate
figure(1)
plot(z,P,'r*');
title("Oxygen Flow Rate Through Fiber Bundle");
xlabel("Membrane Length (m)");
ylabel("O2 Flow Rate (mol/hr)");
figure(2)
plot(z,P1,'m*');
title("Nitrogen Flow Rate Through Fiber Bundle");
xlabel("Membrane Length (m)");
ylabel("N2 Flow Rate (mol/hr)");
figure(3)
plot(z,P2,'b--O');
title("Pressure Drop Across Fiber Bore");
xlabel("Membrane Length (m)");
ylabel("Pressure (bar)");
figure(4)
plot(z,P3,'-O','Color',[0.9290, 0.6940, 0.1250]);
title("Volumetric Flow Rate Through Fiber Bundle");
xlabel("Membrane Length (m)");
ylabel("Volumetric Flow Rate (bar)");


end



