offsets = {'_full','_half'};

for offset = 1:2
    
    leg = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};

    models = {'analysis model','synthesis model'};
    colors = [  1  0  0;...
                0  0  0;...
                0  1  1;...
                0  0  1;...
                0  0 .5];
    widths = [ 1 .5 .5 .5 .5 ];
    gaps = 5:5:50;
    figure
    for model = 1:2

        subplot(1,2,model)
        hold on
        if model == 1
            load(['global_test_analysis',offsets{offset}])
        else
            load(['global_test_synthesis',offsets{offset}])
        end

        for weighting = 1:5
            avg_per_algo = mean(SNRs.weighted(:,:,:,weighting),1,'omitnan');
            avg_per_algo_and_gap = mean(avg_per_algo,3);
            plot(gaps,avg_per_algo_and_gap(:),'color',colors(weighting,:),...
                'linewidth',widths(weighting))
        end

        if model == 1
            load(['global_test_synthesis',offsets{offset}])
        else
            load(['global_test_analysis',offsets{offset}])
        end
        avg_per_algo = mean(SNRs.weighted(:,:,:,1),1,'omitnan');
        avg_per_algo_and_gap = mean(avg_per_algo,3);
        plot(gaps,avg_per_algo_and_gap(:),':','color',colors(1,:),...
            'linewidth',widths(1))

        xlabel('gap length [ms]')
        ylabel('SNR [dB]')
        title(models{model})
        legend([leg,sprintf('%s (%s)',models{3-model},'\texttt{none}')])
        xlim([5 50])
        ylim([0 18])
        set(gca,'xtick',5:5:50)
        set(gca,'ytick',0:2:18)
        grid on
        box on
        if offset == 1
            sgtitle('Results of weighting with full offset')
        else
            sgtitle('Results of weighting with half offset')
        end
    end
end