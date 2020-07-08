%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            plotting gradual                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

fractions = 4:4:20;
gaps = 5:5:50;

colors =  [ 0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5;...
           .5 .5 .5 ];

models  = {'_analysis_','_synthesis_'};
offsets = {'full','half'};
weightingtypes = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};

for offset = 1:2
    figure
    for model = 1:2
        load(['global_test',models{model},offsets{offset}])
        for weighting = 1:5
            subplot(2,5,(model-1)*5 + weighting)
            hold on
            for i = 1:length(fractions)         
                avg_per_algo = mean(SNRs.gradual(:,:,:,weighting,i),1,'omitnan');
                avg_per_algo_and_gap = mean(avg_per_algo,3);
                plot(gaps,avg_per_algo_and_gap(:),'linestyle','-','color',colors(i,:))
            end
            avg_per_algo = mean(SNRs.weighted(:,:,:,weighting),1,'omitnan');
            avg_per_algo_and_gap = mean(avg_per_algo,3);
            plot(gaps,avg_per_algo_and_gap(:),'linestyle','-','color','r')
            xlabel('gap length [ms]')
            ylabel('SNR [dB]')
            if model == 1
                title(['analysis model (',weightingtypes{weighting},')'])
            else
                title(['synthesis model (',weightingtypes{weighting},')'])
            end
            xlim([20 50])
            ylim([0 10])
            set(gca,'xtick',20:5:50)
            set(gca,'ytick',0:2:10)
            grid on
            box on
            legend('$$r=h/4$$','$$r=h/8$$','$$r=h/12$$','$$r=h/16$$','$$r=h/20$$','reference','location','southwest')
        end
    end
    if offset == 1
        sgtitle('Results of gradual inpainting with full offset')
    else
        sgtitle('Results of gradual inpainting with half offset')
    end
end