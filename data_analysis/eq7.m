%make graph of yields for each ground state for each sim type
clear;

%order: flower 


pf = @(x) 0.027 * exp(x) ./ (1+0.027*exp(x));
pturtle = @(x)  0.43 ./ (1+0.027*exp(x));
ptree = @(x) 0.31 ./ (1+0.027*exp(x));
pboat = @(x)  0.26 ./ (1+0.027*exp(x));

E = [1;4;8];
boat = pboat(E);
turtle = pturtle(E);
tree = ptree(E);
flower = pf(E);

%make a bar graph of the data
figure(1)
hold on
bar([flower, turtle, tree, boat]','BarWidth', 1.0);
ax = gca;
ax.FontSize = 16;

%labeling 
set(gca,'XTick',[])
ylabel("Equilibrium Probability")
legend(["E = 1 kT", "E = 4 kT", "E = 8 kT"], 'location','northeast' )