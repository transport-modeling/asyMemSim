function [P] = run_ODE15s(  )
 %ode45 solver
zspan=[0 1]; %meter
W0=[0.042378;0.159423;0.4158;0;5]; %based on concentraion of air at 5 bar, and Q = [L/hr] but # based on 0.01ft/s through 2450-fibre bundle

[z,W_s]=ode15s(@(z,W) mem_gas_sep(z,W),zspan,W0);

P = 2450.*W_s(:,1).*W_s(:,3); %o2
P1 = 2450.*W_s(:,2).*W_s(:,3); %n2
P2 = 0.08314*298.*(W_s(:,1)+W_s(:,2)); %delta P
P3 = 2450*(W_s(:,3)); %volumetric flow rate
P4 = (W_s(:,5)); %volumetric flow rate
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
ylabel("Volumetric Flow Rate (L/hr)");
figure(5)
plot(z,P4,'-O','Color',[0.9290, 0.6940, 0.1250]);
title("Volumetric Flow Rate Through Fiber Bundle");
xlabel("Membrane Length (m)");
ylabel("Volumetric Flow Rate (L/hr)");


end



