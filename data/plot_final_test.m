%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    plotting the average performance                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for offset = 1:2
    
    if offset == 1
        load('final_test_full');
    else
        load('final_test_half');
    end
    algos = {'syn. $$\ell_1$$ (\texttt{none})',...
             'ana. $$\ell_1$$ (\texttt{none})',...
             'syn. $$\ell_1$$ (\texttt{norm})',...
             'ana. $$\ell_1$$ (\texttt{energy})',...
             'syn. $$\ell_1$$ (\texttt{iterative})',...   
             'ana. $$\ell_1$$ (\texttt{iterative})',...
             'ana. $$\ell_1$$ (\texttt{energy}) + tdc',...
             'S-SPAIN H',...
             'A-SPAIN',...
             'OMP',...
             'Janssen'};

    gaps = 5:5:50;     
    figure

    %            DR        CP        wDR       wCP       reDR      reCP      wCPtdc      S-SPAIN   A-SPAIN,  OMP       Janssen
    colors  = [ .5 .5 .5; .5 .5 .5;  0  0  0;  0  0  0;  0  1  1;  0  1  1;  0  0  0;  0  0  1;  0  0  1;  1  0  0;  1  0  0 ];
    styles  = { '--',      '-',      '--',     '-',      '--',     '-',      '-.',     '--',     '-',      ':',     '-'     };
    markers = { 'none',    'none',   'none',   'none',   'none',   'none',   'none',   'none',   'none',   'none',   'none'  };
    widths  = [ 1,         1,        1,        1,        1,        1,        1,        1,        1,        1,        1       ];
    
    % which metrics to plot
    metrics_l = [ 1,... % SNR
                  0,... % fullSNR
                  1,... % PemoQ ODG
                  0,... % PemoQ PSM
                  0,... % PemoQ PSMt
                  0,... % PEAQ ODG
                  ];

    metrics = {'SNR [dB]','fullSNR [dB]','PemoQ ODG','PemoQ PSM','PemoQ PSMt','PEAQ ODG'};
              
    switch sum(metrics_l)
        case 1
            a = 1; b = 1;
        case 2
            a = 1; b = 2;
        case 3
            a = 1; b = 3;
        case 4
            a = 2; b = 2;
        case {5, 6}
            a = 2; b = 3;
    end           
    i = 0; % counter of the plots
    
    %% plotting the selected metrics
    for j = 1:6
        if metrics_l(j)
            i = i + 1;
            subplot(a,b,i)
            hold on
            for algo = 1:length(algos)

                % reading the data, computing mean values and plotting
                switch j
                    case 1
                        avg_per_algo = mean(SNRs(:,algo,:,:),1,'omitnan');
                        avg_per_algo_and_gap = mean(avg_per_algo,4);   
                        plot(gaps,avg_per_algo_and_gap(:),...
                            'color',colors(algo,:),...
                            'linestyle',styles{algo},...
                            'marker',markers{algo},...
                            'linewidth',widths(algo))
                    case 2
                        avg_per_algo = mean(fullSNRs(:,algo,:),1,'omitnan'); 
                    case 3
                        avg_per_algo = mean(PemoQ.ODGs(:,algo,:),1,'omitnan');
                    case 4
                        avg_per_algo = mean(PemoQ.PSMs(:,algo,:),1,'omitnan');
                    case 5
                        avg_per_algo = mean(PemoQ.PSMts(:,algo,:),1,'omitnan');
                    case 6
                        avg_per_algo = mean(PEAQ(:,algo,:),1,'omitnan');
                end                
                if j > 1
                   plot(gaps,avg_per_algo(:),...
                        'color',colors(algo,:),...
                        'linestyle',styles{algo},...
                        'marker',markers{algo},...
                        'linewidth',widths(algo))
                end
                
            end
            
            % legend and axes properties
            legend(algos,'location','southwest','numcolumns',1)
            xlabel('gap length [ms]')
            ylabel(metrics{j})
            xlim([5 50])
            grid on
            box on
            ax = gca;
            ax.XTick = 5:5:50;
            switch j
                case {1, 2}
                    ylim([0 18])
                    ax.YTick = 0:2:18;
                case {3, 6}
                    ylim([-4 0])
                case {4, 5}
                    ylim([0 1])
            end
            
        end
    end

    %% sgtitle
    if offset == 1
        sgtitle('Overall comparison of the methods with full offset');
    else
        sgtitle('Overall comparison of the methods with half offset');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 plotting the performance for each signal                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigs = { 'a08_violin',...
         'a16_clarinet',...
         'a18_bassoon',...
         'a25_harp',...
         'a35_glockenspiel',...
         'a41_celesta',...
         'a42_accordion',...
         'a58_guitar_sarasate',...
         'a60_piano_schubert',...
         'a66_wind_ensemble_stravinsky' };

for offset = 1:2
    
    if offset == 1
        load('final_test_full');
    else
        load('final_test_half');
    end
    
    for metric = find(metrics_l)
        
        switch metric
            case 1
                data = SNRs;
            case 2
                data = fullSNRs;
            case 3
                data = PemoQ.ODGs;
            case 4
                data = PemoQ.PSMs;
            case 5
                data = PemoQ.PSMts;
            case 6
                data = PEAQ;
        end
        
        figure

        for signal = 1:10
            
            subplot(2,5,signal)
            hold on
            
            % reading the data, computing mean values if necessary and plotting
            if metric == 1
                plotdata = mean(data(signal,:,:,:),4,'omitnan');
            else
                plotdata = data(signal,:,:);
            end            
            for algo = 1:length(algos)
                plotvector = plotdata(:,algo,:);
                plot(gaps,plotvector(:),...
                    'color',colors(algo,:),...
                    'linestyle',styles{algo},...
                    'marker',markers{algo},...
                    'linewidth',widths(algo))
            end
            
            % legend and axes properties
            if signal == 10
                legend(algos,'position',[0 0.475 1 0.05],'orientation','horizontal')
            end
            %legend(algos,'location','south')
            xlabel('gap length [ms]')
            ylabel(metrics{metric})
            xlim([5 50])
            grid on
            box on
            ax = gca;
            ax.XTick = 5:5:50;
            switch metric
                case {1, 2}
                    ylim([0 24])
                    ax.YTick = 0:4:24;
                case {3, 6}
                    ylim([-4 0])
                case {4, 5}
                    ylim([0 1])
            end            
            title(sigs{signal},'interpreter','none')
            
        end
        
        %% sgtitle
        if offset == 1
            sgtitle('Comparison of the methods with full offset for each signal');
        else
            sgtitle('Comparison of the methods with half offset for each signal');
        end
    
    end
end