%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      INPAINTING METHODS COMPARISON                      %
%                    testing all the possible variants                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

close all
clear
clc

rng(0)

variant = 'analysis';  % analysis / synthesis
off     = 'full';      % full / half

%% paths
addpath('reweighted l1 relaxation');
load('EBU_SQAM.mat');

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
     
gaps = 5:5:50;
n = 8; % number of gaps

%% fields for the results

% weighted
weighting       = {'none','supp','abs','norm','energy'};
SNRs.weighted   = NaN(length(sigs), length(gaps), n, length(weighting));

% reweighted
SNRs.reweighted = NaN(length(sigs), length(gaps), n);

% gradual
fractions       = 4:4:20;
SNRs.gradual    = NaN(length(sigs), length(gaps), n, length(weighting), length(fractions));

% tdc
addgaps         = 2:2:10;
SNRs.tdc        = NaN(length(sigs), length(gaps), n, length(weighting), length(addgaps));

totalvariants   = length(weighting)*length(fractions) + 1 + length(weighting) + length(weighting)*length(fractions);

load(['data/global_test_',variant,'_',off])

%% main cycles
for signum = 1:length(sigs)
    
for gapnum = 1:length(gaps)

fprintf('\nSignal number: %d/%d',signum,length(sigs));
fprintf('\nGap length number: %d/%d\n\n',gapnum,length(gaps));
    
%% loading signal                                              
signame = sigs{signum};
signal  = eval(signame);

fprintf('Signal: %s\n',signame);
fprintf('Sampling rate: %d Hz\n',fs);

%% transform setting
% window length approximately 64 ms + divisible by 4
w = 2800;
a = w/4;
M = w;

% frame for methods based on convex relaxation
F = frametight(frame('dgt',{'hann',w},a,M,'timeinv'));

%% signal degradation
gap_length  = gaps(gapnum);
h           = round(fs*gap_length/1000); % gap length in samples 
full.length = length(signal);
full.mask   = true(full.length,1);
solution    = signal;

%% segment setting
not_degraded   = 0.5; % length of reliable part at the start / end of the signal (in seconds)
segment.length = round((length(signal)-2*not_degraded*fs)/n);

%% segment processing
soft = @(x,gamma) sign(x) .* max(abs(x)-gamma, 0);
for i = 1:n
    variantcounter = 0;
    fprintf('\nGap number %d of %d.\n\n',i,n);
    idxs = round(not_degraded*fs)+((i-1)*segment.length+1:i*segment.length);
    
    % degrading signal
    s = (w+1) + rand()*(segment.length-2*w-h); % gap start in samples
    
    if ~isnan(SNRs.tdc(signum,gapnum,i,end,end))
        continue
    end
    
    s = round(s);
    f = s + h - 1; % gap end in samples
    segment.mask = true(segment.length,1); % indicator of the reliable samples
    segment.mask(s:f) = 0;
    segment.data = signal(idxs);
    segment.max = max(abs(segment.data));
    segment.data = segment.data/segment.max;
    segment.gapped = segment.data.*segment.mask; % lost samples set to 0 
    full.mask(idxs) = segment.mask;
    
    % shortening the signal
    [q,~,~,~,~,~,~,~,U,V,L] = min_sig_supp_2(w,a,M,s,f,segment.length,1,offset(s,f,a,off));
    origL = L;
    if L < framelength(F,L)
        L = framelength(F,L);
    end
    Q = q + L - 1;
    if Q <= segment.length
        ssegment.mask = segment.mask(q:Q);
        ssegment.gapped = segment.gapped(q:Q);
        ssegment.data = segment.data(q:Q);
    else
        ssegment.mask = segment.mask(q:end);
        ssegment.gapped = segment.gapped(q:end);
        ssegment.data = segment.data(q:end);
        ssegment.mask(end+1:L) = true;
        ssegment.gapped(end+1:L) = 0;
        ssegment.data(end+1:L) = 0;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           gradual inpainting                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for weightcounter = 1:length(weighting)
        weight = weighting{weightcounter};
        fractioncounter = 0;
        for fraction = fractions
            fractioncounter = fractioncounter + 1;
        
            gapindexes = find(~ssegment.mask);
            gradual.s = gapindexes(1);
            gradual.f = gapindexes(end);
            gradual.r = ceil(h/fraction);
            gradual.mask = ssegment.mask;
            gradual.gapped = ssegment.gapped;
            gradual.atoms = syn_atoms(F);

            if strcmp(variant,'analysis')

                % not changing parameters
                CPparam.g = @(x) 0; % it is actually the indicator function... 
                CPparam.K = @(x) frana(F,x);
                CPparam.K_adj = @(x) frsyn(F,x);
                CPparam.dim = length(frsyn(F,frana(F,ssegment.gapped)));

                % the main cycle
                while gradual.s <= gradual.f
                    % computing the weights for l1 relaxation
                    if strcmp(weight,'none')
                        gradual.wts = 1;
                    else
                        gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, weight);
                    end
                        
                    % changing parameters
                    CPparam.f = @(x) norm(gradual.wts.*x,1);
                    CPparam.prox_f = @(x,gamma) soft(x,gradual.wts*gamma);  
                    CPparam.prox_g = @(x,gamma) proj_time(x,gradual.mask,gradual.gapped);

                    % algorithm
                    [ recovered_gradual, ~, ~ ] = ChambollePock(CPparam,[]);
                    recovered_gradual = real(recovered_gradual);

                    % s, f update
                    gradual.s = gradual.s + gradual.r;
                    gradual.f = gradual.f - gradual.r;
                    gradual.mask = true(L,1);
                    gradual.mask(gradual.s:gradual.f) = false;
                    gradual.gapped = gradual.mask.*recovered_gradual;
                end

            else

                % not changing parameters        
                DRparam.g = @(x) 0; % it is actually the indicator function...
                c = frana(F,ssegment.gapped);
                DRparam.dim = length(c);

                % the main cycle
                while gradual.s <= gradual.f
                    % computing the weights for l1 relaxation
                    if strcmp(weight,'none')
                        gradual.wts = 1;
                    else
                        gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, weight);
                    end
                    
                    % changing parameters
                    DRparam.f = @(x) norm(gradual.wts.*x,1);
                    DRparam.prox_f = @(x,gamma) soft(x,gradual.wts*gamma);
                    c = frana(F,gradual.gapped);
                    DRparam.prox_g = @(x,gamma) x - frana(F,gradual.mask.*frsyn(F,x)) + c;

                    % algorithm
                    [ x_hat, ~, ~ ] = DouglasRachford(DRparam,[]);
                    x_hat = DRparam.prox_g(x_hat);
                    recovered_gradual = frsyn(F,x_hat);
                    recovered_gradual = real(recovered_gradual);

                    % s, f update
                    gradual.s = gradual.s + gradual.r;
                    gradual.f = gradual.f - gradual.r;
                    gradual.mask = true(L,1);
                    gradual.mask(gradual.s:gradual.f) = false;
                    gradual.gapped = gradual.mask.*recovered_gradual;
                end

            end
            
            % updating the segment solution
            segment.solution = segment.gapped;
            segment.solution(q:q+origL-1) = recovered_gradual(1:origL);

            % updating the global solution
            solution(idxs) = segment.solution*segment.max;

            % SNR
            SNRs.gradual(signum,gapnum,i,weightcounter,fractioncounter)...
                = snr_n(signal(idxs).*(~segment.mask),solution(idxs).*(~segment.mask));
            variantcounter = variantcounter + 1;
            fprintf('   Done: %2d of %d variants.\n',variantcounter,totalvariants)
            
        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          iterative reweighting                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % parameters of the main cycle
    RWparamsolver.maxit = 10;
    RWparamsolver.epsilon = 1e-3;
    RWparamsolver.delta = 0.01;

    % transform settings
    param.type = variant;
    param.F = F;
    param.offset = off;   
    param.weighting = 'none';
    param.reweighting = true;

    % algorithm
    segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], RWparamsolver);

    % updating the global solution
    solution(idxs) = segment.solution(1:segment.length)*segment.max;
    
    % SNR
    SNRs.reweighted(signum,gapnum,i)...
            = snr_n(signal(idxs).*(~segment.mask),solution(idxs).*(~segment.mask));
    variantcounter = variantcounter + 1;
    fprintf('   Done: %2d of %d variants.\n',variantcounter,totalvariants)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                     direct time-domain compensation                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for weightcounter = 1:length(weighting)
        weight = weighting{weightcounter};
        gapcounter = 0;
        for addgap = addgaps
            gapcounter = gapcounter + 1;
       
            % (1) taking the basic reconstruction of the original gap
        
            % transform settings
            param.type = variant;
            param.F = F;
            param.offset = off;   
            param.weighting = weight;
            param.reweighting = false;

            % algorithm
            segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

            % updating the global solution
            solution(idxs) = segment.solution(1:segment.length)*segment.max;

            % (2) direct time domain compensation (working with the whole signal)

            % ensuring only the single gap
            TDCmask = true(length(solution),1);
            TDCmask(idxs) = segment.mask;

            % TDC parameters
            TDCparam.type = variant;
            TDCparam.F = F;
            TDCparam.offset = off;
            TDCparam.weighting = weight;
            TDCparam.reweighting = false;
            TDCparamsolver.gaps = addgap;
            TDCparamsolver.segs = 10;      
            TDCparamsolver.shift = w/2;
            TDCparamsolver.lens = h/4;

            % compensation
            solution = tdc( solution, TDCmask, TDCparam, [], [], TDCparamsolver );
            
            % SNR
            SNRs.tdc(signum,gapnum,i,weightcounter,gapcounter)...
                = snr_n(signal(idxs).*(~segment.mask),solution(idxs).*(~segment.mask));
            variantcounter = variantcounter + 1;
            fprintf('   Done: %2d of %d variants.\n',variantcounter,totalvariants)
            
        end
        
    end
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                weighting                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for weightcounter = 1:length(weighting)
        weight = weighting{weightcounter};
        
        % settings
        param.type = variant;
        param.F = F;
        param.offset = off;   
        param.weighting = weight;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution(idxs) = segment.solution(1:segment.length)*segment.max;
    
        % SNR        
        SNRs.weighted(signum,gapnum,i,weightcounter)...
            = snr_n(signal(idxs).*(~segment.mask),solution(idxs).*(~segment.mask));
        variantcounter = variantcounter + 1;
        fprintf('   Done: %2d of %d variants.\n',variantcounter,totalvariants)
        
    end
    
    save(['data/global_test_',variant,'_',off],'SNRs')

end % i

end % gapnum

end % signum