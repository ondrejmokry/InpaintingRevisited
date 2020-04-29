weightingtypes = {'none','supp','abs','norm','energy'};
offsets = {'','_half'};

for offset = 1:2
    
    for wts = 1:5
        load(strcat('SNRs_DRCP_',weightingtypes{wts},offsets{offset}));
        eval(strcat(weightingtypes{wts},'.SNRs = SNRs;'));
    end

    legenda = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};

    models = {'synthesis model','analysis model'};
    colors = [  1  0  0;...
                0  0  0;...
                0  1  1;...
                0  0  1;...
                0  0 .5];
    widths = [ 1 .5 .5 .5 .5 ];
    gaps = 5:5:50;
    figure;
    for algo = 1:2

        subplot(1,2,algo);

        for wts = 1:5
            avg_per_algo = mean(eval(strcat(weightingtypes{wts},'.SNRs(:,algo,:,:)')),1);
            avg_per_algo_and_gap = mean(avg_per_algo,4);
            plot(gaps,avg_per_algo_and_gap(:),'color',colors(wts,:),...
                'linewidth',widths(wts));
            hold on;
        end

        if algo == 1
            avg_per_algo = mean(none.SNRs(:,2,:,:),1);
            cislo = 2;
        else
            avg_per_algo = mean(none.SNRs(:,1,:,:),1);
            cislo = 1;
        end
        avg_per_algo_and_gap = mean(avg_per_algo,4);
        plot(gaps,avg_per_algo_and_gap(:),':','color',colors(1,:),...
            'linewidth',widths(1));

        xlabel('gap length [ms]');
        ylabel('SNR [dB]');
        title(models{algo});
        legend([legenda,sprintf('%s (%s)',models{cislo},'\texttt{none}')]);
        ylim([0 18]);
        xlim([5 50]);
        grid on;
        if offset == 1
            sgtitle('Results of weighting with full offset');
        else
            sgtitle('Results of weighting with half offset');
        end
    end

end