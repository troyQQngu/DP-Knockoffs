figure
plot(eps_list,power_list,LineStyle='-',LineWidth=2,Color=[0 0 1])
hold on
plot(eps_list,power_th_list,LineStyle=':',LineWidth=2,Color=[0 0.5 0.75])
plot(eps_list,fdr_list,LineStyle='-',LineWidth=2,Color=[1 0 0])
plot(eps_list,fdr_th_list,LineStyle=':',LineWidth=2,Color=[0.75 0.5 0 ])
ylim([0 1])
xlabel("epsilon")
ylabel("Power/FDR")
legend("Power","Theoretical Power","FDR","Theoretical FDR")
title("Power/FDR tradeoff with Privacy(epsilon)")
