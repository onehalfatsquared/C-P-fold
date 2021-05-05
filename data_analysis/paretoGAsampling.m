%analyze the pareto data from the sampling algorithm
clear; clf;

filename = "paretoGAsampling.txt";
filename = "paretoGAtriangle2harmonic.txt";
P = importdata(filename);
[m,n] = size(P);
x = P(:,1);
y = P(:,2);

scatter(x,y)
xlabel("1-exp(dE_{min})")
ylabel("exp(<E_{misfold}>)")

params1 = P(3,3:end);
same = sum(params1 == P(:,3:end),2);
matching = find(same == n-2);

hold on
scatter(x(matching),y(matching),'b','fill')
axis([0 1 0 1])



cfile = "paretoCorrelations.txt";
P = importdata(cfile);

%split data
eq = P(:,1); rate = P(:,2);
xtest = P(:,3); ytest = P(:,4);
%plot(xtest,ytest)
