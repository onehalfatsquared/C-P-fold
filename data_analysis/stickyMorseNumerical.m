%evaluate the integral for the sticky parameter numerically
clear;

rho = 40;
up = 1.4;

E_disc = linspace(0.001,12,1000);
for i = 1:1000
    E = E_disc(i);
    
    Z = @(r) morse(r,E,rho);
    I = @(r) exp(-Z(r));
    
    kappaNum(i) = integral(I,0,up,'RelTol',0,'AbsTol',1e-12)-(up-1);
    kappaFormula(i) = sqrt(2*pi) * exp(E) / sqrt(2*rho*rho*E);
end


figure(1)
plot(E_disc,log(kappaFormula),E_disc,log(kappaNum))
legend("Laplace Method", "Numerical Integration")
xlabel("Well-Depth")
ylabel("log(Sticky Parameter)")

figure(2)
q = 200;
plot(E_disc(1:q),kappaNum(1:q));














function u = morse(r,E,rho)
    Y = exp(-rho*(r-1));
    u = E*(Y.*Y-2*Y);
end