gap_lengths = 5:5:50;

leg = {'\texttt{none}, half offset',...
       '\texttt{none}, full offset',...
       '\texttt{iterative}, half offset',...
       '\texttt{iterative}, full offset',...
       };

colors = [  0  0  0;...
            0  0  1;...
            0  0  0;...
            0  0  1 ];    
linestyles = {'-','-','--','--'};

figure
for ana = 1:2
    subplot(1,2,ana)
    hold on
    
    for algo = 1:4
        
        switch algo
            case 1 % no weighting, half
                load('results_half')
                data = ana;
            case 2 % no weighting, full
                load('results')
                data = ana;
            case 3 % reweighting, half
                load('results_half')
                data = 2+ana;
            case 4 % reweighting, full
                load('results')
                data = 2+ana;
        end
        
        % SNR measured from single gap
        avg_per_algo = mean(SNRs,1,'omitnan');
        avg_per_algo_and_gap = mean(avg_per_algo,4,'omitnan');
        
        plot(gap_lengths,squeeze(avg_per_algo_and_gap(1,data,:)),'color',colors(algo,:),...
            'linestyle',linestyles{algo})
    end
    legend(leg)
    xlabel('gap length [ms]')
    ylabel('SNR [dB]')
    ylim([0 20])
    xlim([5 50]) 
    grid on
    box on
    if ana == 1
        title('synthesis model')
    else
        title('analysis model')
    end
    sgtitle('Results of iterative reweighting')
end