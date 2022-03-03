function [] = PlotPolyCdfDemonstration3_2D(xTrain,xTildeTrain)
[n, dim] = size(xTrain);
assert(dim == 2);
figure('Name', 'Demonstrate T #3');
subplot(1,2,1)
scatter(xTrain(:,1),xTrain(:,2),50,1:n,'filled');
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
subplot(1,2,2)
scatter(xTildeTrain(:,1),xTildeTrain(:,2),50,1:n,'filled');
xlabel('$\tilde{x}_1$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$\tilde{x}_2$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
end

