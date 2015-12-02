% function pronto_res_plot()


clc;close all;
sub_num = size(PRT.model.output.fold,2);
data_plot = zeros(sub_num,3);
for ii = 1:sub_num
    targets = PRT.model.output.fold(1,ii).targets;
    func_val = PRT.model.output.fold(1,ii).func_val;
    b = PRT.model.output.fold(1,ii).b;
    data_plot(ii,:) = [targets func_val b];
end
figure('Color','w')
x1 = data_plot(1:16,2);
y1 = 1:16;%data_plot(1:16,3);
% plot(x1,y1,'s','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',10);
plot(x1,y1,'s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10,'LineWidth',1.5);
hold on;
x2 = data_plot(17:end,2);
y2 = 17:32;%data_plot(17:end,3);
% plot(x2,y2,'d','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
plot(x2,y2,'^','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10,'LineWidth',1.5)
set(gca,'XLim',[-2.5 2.5],'xtick',-2:1:2,'ytick',[],'LineWidth',1.5,'FontSize',10,'FontWeight','bold')

% title('SVM Result','FontName','Arial','FontSize',14);
xlabel('Function Value','FontName','Arial','FontSize',12,'FontWeight','bold')
ylabel('Subjects','FontName','Arial','FontSize',12,'FontWeight','bold')
legend('Control','Patient')
% add vertical line 
x=[0,0];
y=[0,34];
plot(x,y,'--k','LineWidth',2)

a=1;
export_fig(sprintf('plot%d', a),'-a1','-tiff','-r300');

% plot the ROIC curve

targpos = data_plot(:,1) ==1;
fVals = data_plot(:,2);

% function prt_plot_ROC
[y,idx] = sort(fVals);
targpos = targpos(idx);

fp = cumsum(single(targpos))/sum(single(targpos));
tp = cumsum(single(~targpos))/sum(single(~targpos));
figure('Color','w')

plot(fp,tp,'--ks','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',2);
xlabel('False positives (1-specificity)','FontName','Arial','FontSize',12,'FontWeight','bold')
ylabel('True positives (sensitivity)','FontName','Arial','FontSize',12,'FontWeight','bold')

set(gca,'xtick',0:0.2:1,'ytick',0:0.2:1,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');

a=3;
export_fig(sprintf('plot%d', a),'-a1','-tiff','-r300');