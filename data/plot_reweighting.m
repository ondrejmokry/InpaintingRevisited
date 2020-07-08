gaps = 5:5:50;

leg = {'\texttt{none}, full offset',...
       '\texttt{none}, half offset',...
       '\texttt{iterative}, full offset',...
       '\texttt{iterative}, half offset',...
       };

colors = [  0  0  0;...
            0  0  0;...
            0  0  1;...
            0  0  1 ];    
linestyles = {'-','--','-','--'};

figure
for model = 1:2
    subplot(1,2,model)
    hold on
    if model == 1
        str = 'analysis';
    else
        str = 'synthesis';
    end
    
    for algo = 1:4
        
        switch algo
            case 1 % no weighting, full
                load(['global_test_',str,'_full'])
                avg_per_algo = mean(SNRs.weighted(:,:,:,1),1,'omitnan');
            case 2 % no weighting, half
                load(['global_test_',str,'_half'])
                avg_per_algo = mean(SNRs.weighted(:,:,:,1),1,'omitnan');
            case 3 % iterative reweighting, full
                load(['global_test_',str,'_full'])
                avg_per_algo = mean(SNRs.reweighted(:,:,:),1,'omitnan');
            case 4 % iterative reweighting, half
                load(['global_test_',str,'_half'])
                avg_per_algo = mean(SNRs.reweighted(:,:,:),1,'omitnan');
        end
        
        % SNR measured from single gap
        avg_per_algo_and_gap = mean(avg_per_algo,3,'omitnan');
        
        plot(gaps,squeeze(avg_per_algo_and_gap(:)),'color',colors(algo,:),...
            'linestyle',linestyles{algo})
    end
    legend(leg)
    xlabel('gap length [ms]')
    ylabel('SNR [dB]')
    xlim([5 50])
    ylim([0 18])
    set(gca,'xtick',5:5:50)
    set(gca,'ytick',0:2:18)
    grid on
    box on
    if model == 1
        title('analysis model')
    else
        title('synthesis model')
    end
end
sgtitle('Results of iterative reweighting')