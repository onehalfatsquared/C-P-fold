%make plots that include the scatter data for a given state, along with the
%1 partcle values, and the pareto front

clear; clf;

%get the pareto front data
paretoFile = "paretoOcta001001.txt";
P = importdata(paretoFile);

%do the plotting
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%get the default matlab color order
co = lines(7);

%get data
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
%[AA,BB] = fix_symmetry(AA,BB);

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
        elseif (sticky(i,j) < 0.03)
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


%plot
figure(1)
hold on
plot(EQ, w(:,3), 'color',co(3,:), 'linewidth',3) 
plot(EQ, w(:,2), 'color', co(2,:),'linewidth',3)
plot(EQ, w(:,1), 'color',co(1,:),'linewidth',3)


hold off
xlabel("\pi_S")
%xlabel("Rate")
ylabel("Well Depth (kT)");
%axis([0.2,1,0,15])


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
