figure
subplot(3,1,1)
plot(n_list,log(abs(power_list(1,:)-power_th_list(1,:))./power_th_list(1,:)),LineWidth=2)
hold on 
plot(n_list,zeros(1,length(n_list)),Linestyle=":",LineWidth=2)
xlabel("n")
ylabel("Ratio of Power to Theoretical Power")
%ylim([-0.5 5])
title("Convergence of Power with Sample Size(eps = 0.2)")

subplot(3,1,2)
plot(n_list,log(abs(power_list(2,:)-power_th_list(2,:))./power_th_list(2,:)),LineWidth=2)
hold on 
plot(n_list,zeros(1,length(n_list)),Linestyle=":",LineWidth=2)
xlabel("n")
ylabel("Ratio of Power to Theoretical Power")
title("Convergence of Power with Sample Size(eps = 0.25)")
%ylim([-0.5 5])

subplot(3,1,3)
plot(n_list,log(abs(power_list(3,:)-power_th_list(3,:))./power_th_list(3,:)),LineWidth=2)
hold on 
plot(n_list,zeros(1,length(n_list)),Linestyle=":",LineWidth=2)
xlabel("n")
ylabel("Ratio of Power to Theoretical Power")
%ylim([-0.5 5])
title("Convergence of Power with Sample Size(eps = 0.4)")

figure
subplot(3,1,1)
plot(n_list,log(power_list(1,:)./power_th_list(1,:)),LineWidth=2)
hold on 
plot(n_list,zeros(1,length(n_list)),Linestyle=":",LineWidth=2)
xlabel("n")
ylabel("Ratio of Power to Theoretical Power")
ylim([-0.5 5])
xlim([n_list(1),n_list(end)])
title("Convergence of Power with Sample Size(eps = 0.2)")

subplot(3,1,2)
plot(n_list,log(power_list(2,:)./power_th_list(2,:)),LineWidth=2)
hold on 
plot(n_list,zeros(1,length(n_list)),Linestyle=":",LineWidth=2)
xlabel("n")
ylabel("Ratio of Power to Theoretical Power")
title("Convergence of Power with Sample Size(eps = 0.25)")
ylim([-0.5 5])
xlim([n_list(1),n_list(end)])


subplot(3,1,3)
plot(n_list,log(power_list(3,:)./power_th_list(3,:)),LineWidth=2)
hold on 
plot(n_list,zeros(1,length(n_list)),Linestyle=":",LineWidth=2)
xlabel("n")
ylabel("Ratio of Power to Theoretical Power")
ylim([-0.5 5])
xlim([n_list(1),n_list(end)])
title("Convergence of Power with Sample Size(eps = 0.4)")