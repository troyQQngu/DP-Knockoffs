figure
hold on 
plot(n_list,power_list_jlt,LineWidth=2,color=[0,0,0.9],Marker="o")
plot(n_list,power_list_ag(1,:),LineWidth=2,color=[0.8,0.1,0.8],Marker="o")
plot(n_list,power_list_ag(2,:),LineWidth=2,color=[0.1,0.4,0.5],Marker="o")
plot(n_list,power_list_ag(3,:),LineWidth=2,color=[0.2,0.7,0.2],Marker="o")
title('Power Comparison: JLT vs Laplace Mechanism','FontSize',15)
xlabel('Sample size')
ylabel('Power')
legend('JLT', 'AG1', 'AG2','AG3',location='southeast')
exportgraphics(gcf, 'exp24.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none')


