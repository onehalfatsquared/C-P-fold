%plot the pareto fronts for a state for various different species numbers
%for boat cluster
clear;

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%get the default matlab color order
co = lines(7);
temp = co(4,:); co(4,:) = co(7,:); co(7,:) = temp;

figure(1);
hold on

max_types = 3;
state = "Boat";

%get the unique type pareto data
paretoFileU = "pareto" + state + "Unique.txt";
P = importdata(paretoFileU);
P(:,1) = P(:,1) + 0.02;
%P(:,2) = P(:,2) + 0.02;
plotData(P, co(7,:));
yscale = max(P(:,2))+0.2;

%get the 6 type pareto data
if (max_types >=6)
    paretoFile6 = "pareto" + state + "6.txt";
    P = importdata(paretoFile6);
    plotData(P, co(6,:));
end

%get the 5 type pareto data
if (max_types >=5)
    paretoFile5 = "pareto" + state + "5.txt";
    P = importdata(paretoFile5);
    plotData(P, co(5,:));
end

%get the 4 type pareto data
if (max_types >=4)
    paretoFile4 = "pareto" + state + "4.txt";
    P = importdata(paretoFile4);
    plotData(P, co(4,:));
end

%get the 3 type pareto data
if (max_types >=3)
    paretoFile3 = "pareto" + state + "3.txt";
    P = importdata(paretoFile3);
    plotData(P, co(3,:));
end

%get the 2 type pareto front data
paretoFileA = "pareto" + state + "AAABBBB.txt";
paretoFileB = "pareto" + state + "ABABABB.txt";
PA = importdata(paretoFileA);
PB = importdata(paretoFileB);
PA = combineSets(PA,PB);
PB = combineSets(PB,PA);

plotData(PA, co(2,:));
plotData(PB, co(2,:));

%get the 1 particle data
paretoFile1 = "pareto" + state + "1.txt";
P = importdata(paretoFile1);
plotData(P, co(1,:));

%make labels
xlabel("Equilibrium Probability")
ylabel("\tau^{-1}")
axis([0 1.05 0 yscale])
hold off;

%make parameter plots
params = 1;
if (params == 1)
    
    figure(2)
    hold on
    
    
    plotParameters(PA);
    plotParameters(PB);
    xlabel("Equilibrium Probability")
    ylabel("Well Depth (kT)");
    
    
    
end



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

function [s1] = combineSets(p1, p2)
    %combine two pareto sets
    
    remove = [];
    
    for i=1:length(p1)
        x = p1(i,1); y = p1(i,2);
        for j = 1:length(p2)
            x1 = p2(j,1); y1 = p2(j,2);
            if (x1 > x && y1 > y)
                remove = [remove i];
                break;
            end
        end
    end
    
    p1(remove,:) = [];
    s1 = p1;
    
    
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


function E = invertKappa(k, kappaNum, E_disc)
    
    diff = abs(kappaNum - k);
    ind = find(diff == min(diff));
    E = E_disc(ind);
end


function plotParameters(P)
    %plot parameter data in P
    
    eq = P(:,1);
    [m,n] = size(P);
    n = n-2;

    [EQ,I] = sort(eq);
    
    sticky = P(:,3:end);
    sticky = sticky(I,:);
    
    co = lines(3);
    
    
    %fix jumping due to symmetry
    %[AA,BB] = fix_symmetry(AA,BB);
    
    %load in kappa(E) computed numerically to invert to E
    load("kappaNum.mat");
    load("E_disc.mat");
    
    
    %do the inversion
    w = zeros(m,n); 
    for i = 1:m
        for j = 1:n
            if (sticky(i,j) > 1.42e4)
                w(i,j) = 14;
            elseif (sticky(i,j) < 0.035)
                w(i,j) = 0.1;
            else
                w(i,j) = invertKappa(sticky(i,j), kappaNum, E_disc);
            end
        end
    end
    
    %check for jumps in the Pareto front (EX: parallelogram)
    jump_tol = 0.2; jumpIndex = 0;
    for i=1:m-1
        jump = EQ(i+1) - EQ(i);
        if (jump > jump_tol)
            jumpIndex = i;
            break
        end
    end
    
    if (jumpIndex > 0)
        EQ = [EQ(1:jumpIndex); NaN; EQ(jumpIndex+1:end)];
        w = [w(1:jumpIndex,:); NaN*ones(1,3); w(jumpIndex+1:end,:)];
    end
    
    plot(EQ, w(:,1), 'color', co(1,:),'linewidth',3)
    if (max(EQ)<0.4)
        plot(EQ, w(:,3),'--', 'color',co(3,:), 'linewidth',3)
    else
        plot(EQ, w(:,3), 'color',co(3,:), 'linewidth',3)
    end
    plot(EQ, w(:,2), 'color',co(2,:),'linewidth',3)
    
    
end
        
        
        
        