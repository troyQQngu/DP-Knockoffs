% plot result experiment
figure
plot(mu_list, power_list, linewidth = 2, LineStyle="-", Color=[0,0,1])
hold on
plot(mu_list, power_th_list, LineWidth=2, LineStyle=":", Color=[0,0.5,0.75])
plot(mu_list, fdr_list, Linewidth = 2, LineStyle="-", Color=[1,0,0])
plot(mu_list, fdr_th_list, Linewidth = 2, LineStyle=":", Color=[0.75,0.5,0])
plot(mu_list, fdr_est_th_list, LineWidth= 2, LineStyle="--",Color=[0.75,0,0.5])
plot(mu_list,repmat(q,list_len),Linewidth = 2, LineStyle="-.", Color=[0,0.5,0])
title("Compare Power/FDR with Theoretical Values Against Signal Magnitude")
xlabel('Signal Magnitude(mu)')
xlim([0,mu_list(end)+0.01])
xticks(0.01:0.02:0.23)
ylabel('Power/FDR')
ylim([-0.01,1.01])

legend('Power','Theoretical Power','FDR','Theoretical FDR','Theoretical FDP estimate','Target FDR',Location='northwest')