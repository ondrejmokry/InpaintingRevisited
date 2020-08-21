%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              plotting tdc                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

gaps = 5:5:50;

colors =  [ 0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5;...
           .5 .5 .5 ];

models = {'analysis','synthesis'};
offsets = {'_full','_half'};
weightingtypes = {'\texttt{none}','\texttt{supp}','\texttt{abs}','\texttt{norm}','\texttt{energy}'};

for offset = 1:2
    figure
    for model = 1:2
        load(['global_test_',models{model},offsets{offset}])
        for weighting = 1:5
            subplot(2,5,(model-1)*5 + weighting)
            hold on

            for parameter = 2:2:10
                
                if parameter == 2

                    % basic algorithm
                    avg_per_algo = mean(SNRs.weighted(:,:,:,weighting),1,'omitnan');
                    avg_per_algo_and_gap = mean(avg_per_algo,3);
                    plot(gaps,avg_per_algo_and_gap(:),'r')
                    
                end

                % direct time domain compensation
                avg_per_algo = mean(SNRs.tdc(:,:,:,weighting,parameter/2),1,'omitnan');
                avg_per_algo_and_gap = mean(avg_per_algo,3);
                plot(gaps,avg_per_algo_and_gap(:),'color',colors(parameter/2,:))

            end

            legend('reference',...
                '$$\texttt{gaps} = 2$$',...
                '$$\texttt{gaps} = 4$$',...
                '$$\texttt{gaps} = 6$$',...
                '$$\texttt{gaps} = 8$$',...
                '$$\texttt{gaps} = 10$$')

            title([models{model},' model (',weightingtypes{weighting},')'])
            xlabel('gap length [ms]')
            ylabel('SNR [dB]')
            xlim([5 50])
            ylim([0 18])
            set(gca,'xtick',5:5:50)
            set(gca,'ytick',0:2:18)
            grid on
            box on
        end
    end
    
    if offset == 1
        sgtitle('Results of direct time domain compensation with full offset')
    else
        sgtitle('Results of direct time domain compensation with half offset')
    end
    
end