%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      INPAINTING METHODS COMPARISON                      %
%                            multi-gap version                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;

rng(0);

%% paths and settings
addpath('reweighted l1 relaxation');
addpath('SPAIN');
addpath('Janssen');

PQ = exist('PemoQ','dir');
if PQ
    addpath(genpath('PemoQ'));
end

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

% control, which methods to use in the default settings
%                 DR | CP | reDR | reCP | gradual | tdc | SSPAIN H | ASPAIN | Janssen
turnon = logical([ 1    1      1     1          1     1          1        1         1]);

fprintf('Notes:\n')
fprintf('  - LTFAT needs to be running on your PC (http://ltfat.github.io/)\n')
fprintf('  - some basic settings are performed using the command window\n')
fprintf('  - all the details need to be set inside the code demo.m\n')

fprintf('\nAvailable signals:\n')
for i = 1:length(sigs)
    fprintf('  %2d.  %s\n',i,sigs{i})
end

fprintf('\nAvailable gap lengths:\n')
for i = 1:length(gaps)
    fprintf('  %2d.  %2d ms\n',i,gaps(i))
end

prompt = '\nChoose signal number (1-10): ';
signums = input(prompt);

prompt = '\nChoose gap length number (1-10): ';
gapnums = input(prompt);

prompt = '\nChoose number of gaps per one signal (recommended max. 10): ';
n = input(prompt);

fprintf('\nControl the algorithms to be used:\n')
prompt = '  Use defaults? (0/1): ';
defaults = input(prompt);

if ~defaults
    methods = {...
        'DR         ',...
        'CP         ',...
        'reDR       ',...
        'reCP       ',...
        'gradual    ',...
        'tdc        ',...
        'S-SPAIN H  ',...
        'A-SPAIN    ',...
        'Janssen    '};
    for method = 1:9
        prompt = sprintf('    %s (0/1): ',methods{method});
        turnon(method) = logical(input(prompt));
    end
end    

prompt = '\nPlot spectrograms? (0/1): ';
spectrograms = input(prompt);

sig_counter = 0;
gap_counter = 0;

%% inpainting
for signum = signums
    
sig_counter = sig_counter + 1;
    
for gapnum = gapnums

gap_counter = gap_counter + 1;    
    
fprintf('\nSignal number: %d (%d of %d chosen)',signum,sig_counter,length(signums));
fprintf('\nGap length: %d ms (%d of %d chosen)\n\n',gaps(gapnum),gap_counter,length(gapnums));
    
%% loading signal                                              
signame = sigs{signum};
signal = eval(signame);

fprintf('Signal: %s\n',signame);
fprintf('Sampling rate: %d Hz\n',fs);

%% transform setting
% window length approximately 64 ms + divisible by 4
w = 2800;
a = w/4;
M = w;

F = frametight(frame('dgt',{'hann',w},a,M,'timeinv'));

DFTred = 4;

%% signal degradation
gap_length = gaps(gapnum);
h = round(fs*gap_length/1000); % gap length in samples 
full.length = length(signal);
full.mask = true(full.length,1);

%% segment setting
not_degraded = 0.5; % length of reliable part at the start / end of the signal in seconds
segment.length = round((length(signal)-2*not_degraded*fs)/n);

%% some global parameters
weighting = 'norm';            % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'
off = 'half';                                    % 'full' / 'half' / 'none'
gradual.type = 'analysis';                       % 'analysis' / 'synthesis'
TDCparam.type = 'analysis';                      % 'analysis' / 'synthesis'

%% initializing solutions
solution.DR = signal;
solution.CP = signal;
solution.reDR = signal;
solution.reCP = signal;
solution.gradual = signal;
solution.tdc = signal;
solution.SSPAIN_H = signal;
solution.ASPAIN = signal;
solution.Janssen = signal;

%% segment processing
soft = @(x,gamma) sign(x) .* max(abs(x)-gamma, 0);
for i = 1:n
    fprintf('\nGap number %d of %d.\n',i,n);
    idxs = round(not_degraded*fs)+((i-1)*segment.length+1:i*segment.length);
    
    % degrading signal
    fprintf('Degrading signal.\n');      
    s = (w+1) + rand()*(segment.length-2*w-h); % gap start in samples
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
%                        weighted Douglas-Rachford                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(1)
        fprintf('Starting Douglas-Rachford...\n')

        % computing the weights
        if strcmp(weighting,'none')
            wts = 1;
        else
            wts = weights( F, L, ssegment.mask, U, V, syn_atoms(F), weighting);
        end

        % algorithm settings
        DRparam.f = @(x) norm(wts.*x,1);
        DRparam.g = @(x) 0; % it is actually the indicator function...
        DRparam.prox_f = @(x,gamma) soft(x,wts*gamma);
        c = frana(F,ssegment.gapped);
        DRparam.prox_g = @(x,gamma) x - frana(F,ssegment.mask.*frsyn(F,x)) + c;
        DRparam.dim = length(c);

        % solver settings
        DRparamsolver = [];

        % algorithm
        [ x_hat, ~, ~ ] = DouglasRachford(DRparam,DRparamsolver);
        x_hat = DRparam.prox_g(x_hat);
        recovered_DR = frsyn(F,x_hat);
        recovered_DR = real(recovered_DR(1:origL));

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_DR;

        % updating the global solution
        solution.DR(idxs) = segment.solution*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         weighted Chambolle-Pock                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(2)
        fprintf('Starting Chambolle-Pock...\n');

        % computing the weights
        if strcmp(weighting,'none')
            wts = 1;
        else
            wts = weights( F, L, ssegment.mask, U, V, syn_atoms(F), weighting);
        end

        % algorithm settings
        CPparam.f = @(x) norm(wts.*x,1);
        CPparam.g = @(x) 0; % it is actually the indicator function...
        CPparam.prox_f = @(x,gamma) soft(x,wts*gamma);
        CPparam.prox_g = @(x,gamma) proj_time(x,ssegment.mask,ssegment.gapped);
        CPparam.dim = length(frsyn(F,frana(F,ssegment.gapped)));
        CPparam.K = @(x) frana(F,x);
        CPparam.K_adj = @(x) frsyn(F,x);

        % solver settings
        CPparamsolver = [];

        % algorithm
        [ recovered_CP, ~, ~ ] = ChambollePock(CPparam,CPparamsolver);
        recovered_CP = real(recovered_CP(1:origL));

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_CP;

        % updating the global solution
        solution.CP(idxs) = segment.solution*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       reweighted Douglas-Rachford                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(3)
        fprintf('Starting reweighted Douglas-Rachford...\n');

        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;

        % solver settings
        DRparamsolver = [];
        CPparamsolver = [];

        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'norm';
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, CPparamsolver, DRparamsolver, RWparamsolver);

        % updating the global solution
        solution.reDR(idxs) = segment.solution(1:segment.length)*segment.max;    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        reweighted Chambolle-Pock                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(4)
        fprintf('Starting reweighted Chambolle-Pock...\n');

        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;

        % solver settings
        DRparamsolver = [];
        CPparamsolver = [];

        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'energy';
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, CPparamsolver, DRparamsolver, RWparamsolver);

        % updating the global solution
        solution.reCP(idxs) = segment.solution(1:segment.length)*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                gradual                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(5)
        fprintf('Starting gradual...\n');
        gapindexes = find(~ssegment.mask);
        gradual.s = gapindexes(1);
        gradual.f = gapindexes(end);
        gradual.r = ceil(h/4);
        gradual.mask = ssegment.mask;
        gradual.gapped = ssegment.gapped;
        gradual.atoms = syn_atoms(F);

        if strcmp(gradual.type,'analysis')

            % not changing parameters
            CPparam.g = @(x) 0; % it is actually the indicator function... 
            CPparam.K = @(x) frana(F,x);
            CPparam.K_adj = @(x) frsyn(F,x);
            CPparam.dim = length(frsyn(F,frana(F,ssegment.gapped)));

            % solver settings
            CPparamsolver = [];

            % the main cycle
            while gradual.s <= gradual.f
                % computing the weights for l1 relaxation
                gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, weighting);

                % changing parameters
                CPparam.f = @(x) norm(gradual.wts.*x,1);
                CPparam.prox_f = @(x,gamma) soft(x,gradual.wts*gamma);  
                CPparam.prox_g = @(x,gamma) proj_time(x,gradual.mask,gradual.gapped);

                % algorithm
                [ recovered_gradual, ~, ~ ] = ChambollePock(CPparam,CPparamsolver);
                recovered_gradual = real(recovered_gradual);

                % s, f update
                gradual.s = gradual.s + gradual.r;
                gradual.f = gradual.f - gradual.r;
                gradual.mask = true(L,1);
                gradual.mask(gradual.s:gradual.f) = false;
                gradual.gapped = gradual.mask.*recovered_gradual;
            end

            % updating the segment solution
            segment.solution = segment.gapped;
            segment.solution(q:q+origL-1) = recovered_gradual(1:origL);

            % updating the global solution
            solution.gradual(idxs) = segment.solution*segment.max;

        else        

            % not changing parameters        
            DRparam.g = @(x) 0; % it is actually the indicator function...
            c = frana(F,ssegment.gapped);
            DRparam.dim = length(c);

            % solver settings
            DRparamsolver = [];

            % the main cycle
            while gradual.s <= gradual.f
                % computing the weights for l1 relaxation
                gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, weighting);

                % changing parameters
                DRparam.f = @(x) norm(gradual.wts.*x,1);
                DRparam.prox_f = @(x,gamma) soft(x,gradual.wts*gamma);
                c = frana(F,gradual.gapped);
                DRparam.prox_g = @(x,gamma) x - frana(F,gradual.mask.*frsyn(F,x)) + c;

                % algorithm
                [ x_hat, ~, ~ ] = DouglasRachford(DRparam,DRparamsolver);
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

            % updating the segment solution
            segment.solution = segment.gapped;
            segment.solution(q:q+origL-1) = recovered_gradual(1:origL);

            % updating the global solution
            solution.gradual(idxs) = segment.solution*segment.max;

        end     
    end      
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 direct time-domain compensation (tdc)                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    if turnon(6)      
        fprintf('Starting direct time-domain compensation (tdc)...\n');

        % (1) taking the basic reconstruction of the original gap (check weighting!)

        if strcmp(TDCparam.type,'synthesis')
            solution.tdc(idxs) = solution.reDR(idxs);
        else
            solution.tdc(idxs) = solution.reCP(idxs);
        end

        % (2) direct time domain compensation (working with the whole signal)

        % ensuring only the single gap
        TDCmask = true(length(solution.tdc),1);
        TDCmask(idxs) = segment.mask;

        % TDC parameters
        TDCparam.F = F;
        TDCparam.offset = off;
        TDCparam.weighting = 'energy';
        TDCparam.reweighting = false;
        TDCparamsolver.gaps = 4;
        TDCparamsolver.segs = 10;      
        TDCparamsolver.shift = w/2;
        TDCparamsolver.lens = h/4;

        % compensation
        solution.tdc = tdc( solution.tdc, TDCmask, TDCparam, [], [], TDCparamsolver );
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               S-SPAIN H                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(7)
        fprintf('Starting S-SPAIN H...\n');

        SPAINparam.algorithm = 'sspain';
        SPAINparam.w = w;
        SPAINparam.a = a;
        SPAINparam.wtype = 'hann';
        SPAINparam.M = M;
        SPAINparam.Ls = L;
        SPAINparam.mask = ssegment.mask;
        SPAINparam.F = frame('dft');
        SPAINparam.F.redundancy = DFTred;
        SPAINparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPAINparam.F.redundancy-1),1)]);
        SPAINparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPAINparam.F.redundancy);

        % solver settings
        SPAINparamsolver.s = 1; % increment of k
        SPAINparamsolver.r = 1; % every r-th iteration increment k by s   
        SPAINparamsolver.epsilon = 0.1; % stopping criterion of termination function
        SPAINparamsolver.maxit = ceil(floor(SPAINparam.w*SPAINparam.F.redundancy/2+1)*SPAINparamsolver.r/SPAINparamsolver.s); % maximum number of iterations
        SPAINparamsolver.store_snr = 0; 
        SPAINparamsolver.store_obj = 0;
        SPAINparamsolver.f_update = 'H';

        % algorithm    
        [recovered_SSPAIN, ~, ~] = spain_segmentation(ssegment.gapped,SPAINparam,SPAINparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_SSPAIN(1:origL);

        % updating the global solution
        solution.SSPAIN_H(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                A-SPAIN                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(8)
        fprintf('Starting A-SPAIN...\n');

        SPAINparam.algorithm = 'aspain';
        SPAINparam.w = w;
        SPAINparam.a = a;
        SPAINparam.wtype = 'hann';
        SPAINparam.M = M;
        SPAINparam.Ls = L;
        SPAINparam.mask = ssegment.mask;
        SPAINparam.F = frame('dft');
        SPAINparam.F.redundancy = DFTred;
        SPAINparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPAINparam.F.redundancy-1),1)]);
        SPAINparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPAINparam.F.redundancy);

        % solver settings
        SPAINparamsolver.s = 1; % increment of k
        SPAINparamsolver.r = 1; % every r-th iteration increment k by s   
        SPAINparamsolver.epsilon = 0.1; % stopping criterion of termination function
        SPAINparamsolver.maxit = ceil(floor(SPAINparam.w*SPAINparam.F.redundancy/2+1)*SPAINparamsolver.r/SPAINparamsolver.s); % maximum number of iterations
        SPAINparamsolver.store_snr = 0; 
        SPAINparamsolver.store_obj = 0;

        % algorithm    
        [recovered_ASPAIN, ~, ~] = spain_segmentation(ssegment.gapped,SPAINparam,SPAINparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_ASPAIN(1:origL);

        % updating the global solution
        solution.ASPAIN(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                Janssen                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(9)
        fprintf('Starting Janssen...\n');

        Janssenparam.w = w;
        Janssenparam.a = a;
        Janssenparam.wtype = 'hann';
        Janssenparam.Ls = L;
        Janssenparam.mask = ssegment.mask;

        % solver settings
        Janssenparamsolver.Nit = 50; % number of iterations

        % algorithm    
        recovered_Janssen = janssen(ssegment.gapped,Janssenparam,Janssenparamsolver);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_Janssen(1:origL);

        % updating the global solution
        solution.Janssen(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
end

%% figures
fulltime = (1:length(signal))/fs;

% full time
leg = {'original','DR','CP','reDR','reCP','gradual','tdc','S-SPAIN H','A-SPAIN','Janssen'};
figure()
hold on
plot(fulltime,signal)
if turnon(1)    
    plot(fulltime,solution.DR)
end
if turnon(2)
    plot(fulltime,solution.CP)
end
if turnon(3)    
    plot(fulltime,solution.reDR)
end
if turnon(4)
    plot(fulltime,solution.reCP)
end
if turnon(5)
    plot(fulltime,solution.gradual)
end
if turnon(6)
    plot(fulltime,solution.tdc)
end
if turnon(7)
    plot(fulltime,solution.SSPAIN_H)
end
if turnon(8)
    plot(fulltime,solution.ASPAIN)
end
if turnon(9)
    plot(fulltime,solution.Janssen)
end
legend(leg([true turnon]))
xlabel('time [s]')
title(sprintf('signal: %s, gap length: %d ms, full signal',signame,gap_length),'interpreter','none')

% only the gaps
figure()
hold on
plot(signal(~full.mask))
if turnon(1)   
    plot(solution.DR(~full.mask))
end
if turnon(2)
    plot(solution.CP(~full.mask))
end
if turnon(3)   
    plot(solution.reDR(~full.mask))
end
if turnon(4)
    plot(solution.reCP(~full.mask))
end
if turnon(5)
    plot(solution.gradual(~full.mask))
end
if turnon(6)
    plot(solution.tdc(~full.mask))
end
if turnon(7)
    plot(solution.SSPAIN_H(~full.mask))
end
if turnon(8)
    plot(solution.ASPAIN(~full.mask))
end
if turnon(9)
    plot(solution.Janssen(~full.mask))
end
legend(leg([true turnon]))
title(sprintf('signal: %s, gap length: %d ms, only the gaps',signame,gap_length),'interpreter','none')

%% SNRs and ODGs
fprintf('\nSNRs and ODGs computed from all the gaps at once:\n')
if turnon(1)
    fprintf('   DR:          SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.DR(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.DR, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(2)
    fprintf('   CP:          SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.CP(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.CP, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(3)
    fprintf('   reDR:        SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.reDR(~full.mask)))    
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.reDR, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(4)
    fprintf('   reCP:        SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.reCP(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.reCP, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(5)
    fprintf('   gradual:     SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.gradual(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.gradual, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(6)
    fprintf('   tdc:         SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.tdc(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.tdc, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(7)
    fprintf('   S-SPAIN H:   SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.SSPAIN_H(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.SSPAIN_H, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(8)
    fprintf('   A-SPAIN:     SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.ASPAIN(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.ASPAIN, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end
if turnon(9)
    fprintf('   Janssen:     SNR: %5.2f dB\n',snr_n(signal(~full.mask),solution.Janssen(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, solution.Janssen, fs);
        fprintf(repmat('\b', 1, 22))
        fprintf('                ODG: %5.2f\n',ODG)
    end
end

%% spectrograms
if spectrograms
    number = sum(turnon)+1;
    facts = factor(number);
    A = facts(1);
    B = max(A,number/A);
    A = number/B;
    titles = leg([true turnon]);
    fields = fieldnames(solution);
    fields = fields(turnon);
    figure
    subplot(A,B,1)
    sg(signal)
    title(titles{1})
    for spec = 2:number
        subplot(A,B,spec)
        sg(solution.(fields{spec-1}))
        title(titles{spec})
    end
end

end % gapnum

end % signum