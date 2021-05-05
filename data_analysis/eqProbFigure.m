%make bar graph for the eq prob figure
clear;

set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultAxesFontSize',18)

%set the probs
%order - 1 particle, unique particles, ABABAB binary, ABABAB free
parallel = [49, 0.5, 0.5, 0.5]; 
chevron = [40, 0.5, 80, 7];
triangle = [11, 100, 20, 93];

%make a bar graph of the data
figure(1)
bar([parallel; chevron; triangle])

%labeling 
set(gca,'XTick',[])
ylabel("Equilibrium Probability (%)")








%for the 3d system
octa = [5.3, 99.99, 6.6, 90.5]/100;
poly = [94.7, 0.5, 93.3, 7.3]/100; 

figure(2)

%make a bar graph of the data
bar([poly; octa])

%labeling 
set(gca,'XTick',[])
ylabel("Equilibrium Probability")



