%make plots that include the scatter data for a given state, along with the
%1 partcle values, and the pareto front

clear;

params = 1;

%get the pareto front data
paretoFile = "paretoTriangle2.txt";
state = "chevron";
P = importdata(paretoFile);

%get unique data
% uniqueFile = "paretoGA.txt";
% Pu = importdata(uniqueFile);
% Pu(:,1) = Pu(:,1) + 0.02;
% %Pu(:,2) = Pu(:,2) + 0.01;
% 
% %get 3 particle data
% file3 = "paretoTriangle3.txt";
% P3 = importdata(file3);

%get the scatter plot data - 2 types
scatter3file = state + "Scatter3.txt";
scatter3data = importdata(scatter3file); 
prob3 = scatter3data(:,1); rate3 = scatter3data(:,2);

%get the scatter plot data - 1 types
scatter1file =  state + "Scatter1.txt";
scatter1data = importdata(scatter1file); 
prob1 = scatter1data(:,1); rate1 = scatter1data(:,2);










%do the plotting

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%get the default matlab color order
co = lines(7);

%make the figure
figure(1);
hold on;

%plot the scatter 3 first
scatter(prob3, rate3, [], co(1,:));

%plot the scatter 1 next
plot(prob1, rate1, '.-', 'color',co(3,:),'markersize',20,'linewidth',2)
%plot(prob1, rate1,'o', 'color',co(1,:))

%plot unique front
%plotData(Pu,co(4,:));

%plot 3 front
%plotData(P3, co(5,:));

%plot the pareto front
plotData(P,co(2,:));

%make labels
xlabel("Equilibrium Probability")
ylabel("\tau^{-1}")
%axis tight
%axis([0,1,0,inf])


hold off

if (params == 1)
    
    eq = P(:,1);
    rate = P(:,2);
    [m,n] = size(P);
    ind = 1:m;

    [EQ,I] = sort(eq);
    
    aa = P(:,3);
    ab = P(:,4);
    bb = P(:,5);
    
    AA = aa(I);
    AB = ab(I);
    BB = bb(I);
    
    %fix jumping due to symmetry
    [AA,BB] = fix_symmetry(AA,BB);
    
    %load in kappa(E) computed numerically to invert to E
    load("kappaNum.mat");
    load("E_disc.mat");
    
    
    %do the inversion
    sticky = [AA AB BB];
    w = zeros(m,3); 
    for i = 1:m
        for j = 1:3
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
    
    %fix syymmetry issues
    %[w(:,1),w(:,3)] = fix_symmetry(w(:,1),w(:,3));
        
    
    %plot
    
    figure(2)
    hold on
    plot(EQ, w(:,3), 'color',co(1,:), 'linewidth',2)
    plot(EQ, w(:,1), 'color', co(3,:),'linewidth',2)
    plot(EQ, w(:,2), 'color',co(2,:),'linewidth',2)
    
    
    hold off
    xlabel("Equilibrium Probability")
    ylabel("Well Depth (kT)");
    %axis([0.2,1,0,15])
    
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

function plotData(P,color)
    
    [m,n] = size(P);

    eq = P(:,1);
    rate = P(:,2);

    [EQ,I] = sort(eq);
    RATE = rate(I);
    ind = 1:m;
    
    
    
    plot(EQ(ind),RATE(ind),'.', 'color',color,'markersize',20)
end

function E = invertKappa(k, kappaNum, E_disc)
    
    diff = abs(kappaNum - k);
    ind = find(diff == min(diff));
    E = E_disc(ind);
end
