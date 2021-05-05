%make a figure corresponding to the pareto movie for the triangle
clear;

%get the full data set. set parameters related to indexing
P = importdata("triangleMovie2.txt");
[m,n] = size(P);
pop_size = 1000;
gens = 40;

eq = zeros(gens,pop_size);
r = zeros(gens, pop_size);

for i = 1:gens
    ind = ((i-1)*pop_size + 1) : i*pop_size;
    eq(i,:) = P(ind,1);
    r(i,:) = P(ind,2);  
end

genp = 7;
co = jet(genp);

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%first figure, all points, every 2nd generation, color gradient
figure(1)
hold on
for i = 1:genp
    scatter(eq(4*i-3,:),r(4*i-3,:),[],co(i,:));
end
hold off
xlabel("Equilibrium Probability")
ylabel("\tau^{-1}")

%second figure - use less points, show the paths traced out
figure(2)
hold on
for i=1:100:500
    plot(eq(:,i),r(:,i))
end
hold off
xlabel("Equilibrium Probability")
ylabel("\tau^{-1}")

figure(3)
test = 25;
scatter(eq(test,:),r(test,:))
xlabel("Equilibrium Probability")
ylabel("\tau^{-1}")
title("Generation " + num2str(test))
axis([0,1,0,0.45])






