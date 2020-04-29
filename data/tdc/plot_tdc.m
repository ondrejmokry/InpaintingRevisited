gaps = 5:5:50;

colors =  [ 0  0  0;...
            0  1  1;...
            0  0  1;...
            0  0 .5;...
           .5 .5 .5 ];

models = {'ana','syn'};
modelstrings = {'analysis','synthesis'};
offsets = {'','_half'};

for offset = 1:2

    figure
    for model = 1:2
        if model == 1
            weightingtypes = {'none','energy'};
            loc = 'southwest';
            order = [2 4];
        else
            if offset == 1
                weightingtypes = {'none','abs'};
            else
                weightingtypes = {'none','norm'};
            end
            loc = 'northeast';
            order = [1 3];
        end
        for wts = 1:2

            subplot(2,2,order(wts))
            hold on

            for parametr = 2:2:10

                % load the data
                load(sprintf('SNRs_tdc_%s_%s_%d%s',models{model},weightingtypes{wts},parametr,offsets{offset}));

                if parametr == 2

                    % basic algorithm
                    avg_per_algo = mean(SNRs(:,1,:,:),1,'omitnan');
                    avg_per_algo_and_gap = mean(avg_per_algo,4);
                    plot(gaps,avg_per_algo_and_gap(:),'r');
                    
                end

                % direct time domain compensation
                avg_per_algo = mean(SNRs(:,2,:,:),1,'omitnan');
                avg_per_algo_and_gap = mean(avg_per_algo,4);
                plot(gaps,avg_per_algo_and_gap(:),'color',colors(parametr/2,:))

            end

            legend('reference',...
                '$$\texttt{gaps} = 2$$',...
                '$$\texttt{gaps} = 4$$',...
                '$$\texttt{gaps} = 6$$',...
                '$$\texttt{gaps} = 8$$',...
                '$$\texttt{gaps} = 10$$',...
                'Location',loc)

            title(sprintf('model: %s, weighting: \\texttt{%s}',modelstrings{model},weightingtypes{wts}));
            xlabel('gap length [ms]')
            ylabel('SNR [dB]')
            ax = gca;
            ax.XTick = 5:5:50;
            ax.YTick = 0:2:18;

            grid on
            box on
            xlim([5 50])
            ylim([0 18])

        end
    end
    
    if offset == 1
        sgtitle('Results of direct time domain compensation with full offset');
    else
        sgtitle('Results of direct time domain compensation with half offset');
    end
    
end