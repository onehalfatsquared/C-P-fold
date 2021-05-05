%plot the pareto fronts for a state for various different species numbers
%for 6 particle clusters

clear;clf

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%load in kappa(E) computed numerically to invert to E
load("kappaNum.mat");
load("E_disc.mat");

%get the default matlab color order
co = lines(4);

%set the state
state = "Chevron";

%set the interaction of interest
ixn = 3;

%init figures
figure(1);
hold on
figure(2);
hold on


%plot
%get the standard pareto front data
paretoFile1 = "pareto" + state + "ABABAB.txt";
P = importdata(paretoFile1);
[P(:,3),P(:,5)] = fix_symmetry(P(:,3),P(:,5));
plotData(P, co(1,:));
plotParam(P, co(1,:), ixn, kappaNum, E_disc)
yscale(1) = max(P(:,2))+0.1;

%get the perturb 20  pareto front data
paretoFile2 = "perturb" + state + "20.txt";
P = importdata(paretoFile2);
[P(:,3),P(:,5)] = fix_symmetry(P(:,3),P(:,5));
plotData(P, co(2,:));
plotParam(P, co(2,:), ixn, kappaNum, E_disc)
yscale(2) = max(P(:,2))+0.1;

%get the perturb 20  pareto front data
paretoFile3 = "perturb" + state + "C.txt";
P = importdata(paretoFile3);
[P(:,3),P(:,5)] = fix_symmetry(P(:,3),P(:,5));
plotData(P, co(3,:));
plotParam(P, co(3,:), ixn, kappaNum, E_disc)
yscale(3) = max(P(:,2))+0.1;


%make labels for figure 1
figure(1)
xlabel("\pi_S")
ylabel("\tau^{-1}")
axis([0 1.0 0 max(yscale)])
legend("Unperturbed", "Rate Perturbed", "Probability Perturbed")



%make labels for figure 2
figure(2)
xlabel("\pi_S")
ylabel("Well Depth (kT)");
axis([0,1,0,15])
legend("Unperturbed", "Rate Perturbed", "Probability Perturbed")




function plotData(P,color)
    %plot data in P
    
    eq = P(:,1);
    rate = P(:,2);

    [EQ,I] = sort(eq);
    RATE = rate(I);
    [EQ,RATE] = getPareto(EQ,RATE);
    
    figure(1);
    plot(EQ,RATE,'.', 'color',color,'markersize',20)
end

function plotParam(P,color,ixn,kappaNum, E_disc)
    %plot data in P
    
    eq = P(:,1);
    sticky = P(:,ixn+2);

    [EQ,I] = sort(eq);
    STICKY = sticky(I);
    E = getE(STICKY, kappaNum, E_disc);
    
    figure(2);
    plot(EQ,E, 'color',color,'linewidth',3)
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

function E = invertKappa(k, kappaNum, E_disc)
    %invert a kappa value using the data
    diff = abs(kappaNum - k);
    ind = find(diff == min(diff));
    E = E_disc(ind);
end


function w = getE(sticky, kappaNum, E_disc)
    %invert vector of sticky parameters
    
    m = length(sticky);
    w = zeros(1,m);
    
    for i = 1:m
        if (sticky(i) > 1.42e4)
            w(i) = 14;
        elseif (sticky(i) < 0.03)
            w(i) = 0.1;
        else
            w(i) = invertKappa(sticky(i), kappaNum, E_disc);
        end
    end
    
end

function [x,y] = fix_symmetry(a,b)
    %if two curves swap b/c of a particle syymetry, make them continuous
    
    n = length(a);
    x = zeros(n,1);
    y = zeros(n,1);
    
    for i=1:n
        x(i) = min(a(i), b(i));
        y(i) = max(a(i), b(i));
    end
end
    
    
    
    


