fractions = 4:4:20;
gaps = 5:5:50;

styles = {'--','-'};
colors =  [ 0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5;...
           .5 .5 .5 ];

strings = {'_','_syn_'};

for offset = 1:2
    figure
    for model = 1:2
        subplot(1,2,3-model)
        hold on
        for i = 1:length(fractions)
            if offset == 1
                eval(sprintf('load(''SNRs_gradual%s%d.mat'');',strings{model},fractions(i)));
                avg_per_algo = mean(SNRs(:,2,:,:),1,'omitnan');
            else
                eval(sprintf('load(''SNRs_gradual%s%d_half.mat'');',strings{model},fractions(i)));
                avg_per_algo = mean(SNRs(:,1,:,:),1,'omitnan');
            end            
            avg_per_algo_and_gap = mean(avg_per_algo,4);
            plot(gaps,avg_per_algo_and_gap(:),'linestyle','-','color',colors(i,:));        
        end
        if offset == 1
            if model == 1
                load('../weighting/SNRs_DRCP_energy')
                avg_per_algo = mean(SNRs(:,2,:,:),1,'omitnan');
            else
                load('../weighting/SNRs_DRCP_abs')
                avg_per_algo = mean(SNRs(:,1,:,:),1,'omitnan');
            end
        else
            if model == 1
                load('../weighting/SNRs_DRCP_energy_half')
                avg_per_algo = mean(SNRs(:,2,:,:),1,'omitnan');
            else
                load('../weighting/SNRs_DRCP_norm_half')
                avg_per_algo = mean(SNRs(:,1,:,:),1,'omitnan');
            end
        end        
        avg_per_algo_and_gap = mean(avg_per_algo,4);
        plot(gaps,avg_per_algo_and_gap(:),'linestyle','-','color','r')
        xlabel('gap length [ms]')
        ylabel('SNR [dB]')
        if model == 1
            title('analysis model')
        else
            title('synthesis model')
        end
        grid on
        box on
        xlim([20 50])
        ylim([2 10])
        legend('$$r=h/4$$','$$r=h/8$$','$$r=h/12$$','$$r=h/16$$','$$r=h/20$$','reference','location','southwest')
        if offset == 1
            sgtitle('Results of gradual inpainting with full offset')
        else
            sgtitle('Results of gradual inpainting with half offset')
        end
    end
end