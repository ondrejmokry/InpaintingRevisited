colors = [  1  0  0;...
            0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5];
        
alpha = 0.5;

wtsleg = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};
wts = {'none','supp','abs','norm','energy'};
figure
for model = 1:2
    subplot(1,2,model)
    hold on
    for i = 1:5        
        load(['SNRs_DRCP_',wts{i}])
        x = SNRs(:,model,:,:);
        load(['SNRs_DRCP_',wts{i},'_half'])
        y = SNRs(:,model,:,:);
        h(i) = scatter(x(:),y(:),[],colors(i,:),'MarkerEdgeAlpha',alpha);
    end
    xlabel('SNR [dB] for full offset')
    ylabel('SNR [dB] for half offset')
    grid on
    box on
    if model == 1
        title('synthesis model')
    else
        title('analysis model')
    end
    X = get(gca,'XLim');
    Y = get(gca,'YLim');
    maximum = max([X(:); Y(:)]);
    line([-5, maximum+5],[-5, maximum+5],'Color',[.5 .5 .5],'LineWidth',2)
    xlim([0, maximum])
    ylim([0, maximum])
    legend(h,wtsleg,'location','southeast')
end
sgtitle('Comparison of the two offset approaches')