%plot the pareto fronts for a state for various different species numbers
%for 6 particle clusters
clear;

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%get the default matlab color order
co = lines(4);

figure(1);
hold on

max_types = 3;
state = "Triangle";

%get the unique type pareto data
paretoFileU = "pareto" + state + "Unique.txt";
P = importdata(paretoFileU);
P(:,1) = P(:,1) + 0.02;
%P(:,2) = P(:,2) + 0.02;
plotData(P, co(4,:));
yscale = max(P(:,2))+0.2;

%get the 3 type pareto data
if (max_types >=3)
    paretoFile3 = "pareto" + state + "3.txt";
    P = importdata(paretoFile3);
    plotData(P, co(3,:));
end

%get the 2 type pareto front data
paretoFile2 = "pareto" + state + "2.txt";
P = importdata(paretoFile2);
plotData(P, co(2,:));

%get the 1 particle data
paretoFile1 = "pareto" + state + "1.txt";
P = importdata(paretoFile1);
plotData(P, co(1,:));





%make labels
xlabel("Equilibrium Probability")
ylabel("\tau^{-1}")

axis([0 1.05 0 yscale])


function plotData(P,color)
    %plot data in P
    
    eq = P(:,1);
    rate = P(:,2);

    [EQ,I] = sort(eq);
    RATE = rate(I);
    [EQ,RATE] = getPareto(EQ,RATE);
    
    plot(EQ,RATE,'.', 'color',color,'markersize',20)
end

function [x,y] = getPareto(eq, rate)
    %remove the increasing part of 1 particle data
    
    M = max(rate);
    M_ind = find(rate == M);
    L = length(rate);
    pareto_ind = M_ind:L;
    
    x = eq(pareto_ind);
    y = rate(pareto_ind);
end