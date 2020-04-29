weightingtypes = {'none','supp','abs','norm','energy'};
leg = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};
colors = [  1  0  0;...
            0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5];
widths = [ 1 .5 .5 .5 .5 ];

% load data (gap length 35 ms)
load('weights.mat');

% plot
figure;
for i = 1:5
    subplot(1,5,i);
    ax = plot(eval(weightingtypes{i}),'color',colors(i,:),'linewidth',1.5); hold on
    legend(leg{i})
    xlabel('index atomu $$n$$')
    ylabel('v\''{a}ha $$w_n$$')
    xlim([0 30]*1000)
    ylim([0 1])
    grid on
    grid minor
end
sgtitle('Illustration of the weighting variants')