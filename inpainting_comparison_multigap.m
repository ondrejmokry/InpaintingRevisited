%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      INPAINTING METHODS COMPARISON                      %
%                            multi-gap version                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% algos
% DR ........ (weighted) synthesis l1
% CP ........ (weighted) analysis l1
% reDR ...... reweighted synthesis l1
% reCP ...... reweighted synthesis l1
% gradual ... gradual approach
% tdc ....... direct time domain compensation
% SSPAIN_H .. synthesis variant of SPAIN 
% ASPAIN .... analysis variant of SPAIN
% Janssen ... AR-based algorithm by Janssen et al.

close all;
clear;
clc;

rng(0);

%% paths
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
n = 8; % number of gaps

%              signals       algos  gap length    gap position
SNRs     = NaN(length(sigs), 9,     length(gaps), n);
fullSNRs = NaN(length(sigs), 9,     length(gaps));
ODGs     = NaN(length(sigs), 9,     length(gaps));
PSMs     = NaN(length(sigs), 9,     length(gaps));
PSMts    = NaN(length(sigs), 9,     length(gaps));

for signum = 1:length(sigs)

for gapnum = 1:length(gaps)
        
fprintf('\nSignal number: %d/%d',signum,length(sigs));
fprintf('\nGap length number: %d/%d\n\n',gapnum,length(gaps));
    
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
not_degraded = 0.5; % length of reliable part at the start / end of the signal
segment.length = round((length(signal)-2*not_degraded*fs)/n);

%% some global parameters
weighting = 'none';            % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'    
off = 'half';                                    % 'full' / 'half' / 'none'
gradual.type  = 'analysis';                      % 'analysis' / 'synthesis'
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         weighted Chambolle-Pock                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       reweighted Douglas-Rachford                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        reweighted Chambolle-Pock                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                gradual                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 direct time-domain compensation (tdc)                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               S-SPAIN H                                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                A-SPAIN                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                Janssen                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    %% saving SNRs

    % weighted l1
    SNRs(signum,1,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.DR(idxs).*(~segment.mask));
    SNRs(signum,2,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.CP(idxs).*(~segment.mask));

    % iteratively reweighted l1
    SNRs(signum,3,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.reDR(idxs).*(~segment.mask));
    SNRs(signum,4,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.reCP(idxs).*(~segment.mask));

    % gradual
    SNRs(signum,5,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.gradual(idxs).*(~segment.mask));
    
    % tdc
    SNRs(signum,6,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.tdc(idxs).*(~segment.mask));
    
    % SPAIN
    SNRs(signum,7,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.SSPAIN_H(idxs).*(~segment.mask));
    SNRs(signum,8,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.ASPAIN(idxs).*(~segment.mask));
    
    % Janssen
    SNRs(signum,9,gapnum,i) = snr_n(signal(idxs).*(~segment.mask),solution.Janssen(idxs).*(~segment.mask));
    
end

%% saving SNRs and perceptual similarity measures
measures = {'fullSNRs','ODGs','PSMs','PSMts'};
algos = {'DR','CP','reDR','reCP','gradual','tdc','SSPAIN_H','ASPAIN','Janssen'};
 
for algo = 1:9
    eval(sprintf('fullSNRs(signum,algo,gapnum) = snr_n(signal(~full.mask),solution.%s(~full.mask));',algos{algo}));
    
    fprintf('\nAlgorithm: %s\n',algos{algo})
    fprintf('  SNR:  %5.2f dB\n',fullSNRs(signum,algo,gapnum))
    
    if PQ
        eval(sprintf('[PSM, PSMt, ODG, ~] = audioqual(signal, solution.%s, fs);',algos{algo}));
        fprintf(repmat('\b', 1, 22))

        PSMs (signum,algo,gapnum) = PSM;
        PSMts(signum,algo,gapnum) = PSMt;
        ODGs (signum,algo,gapnum) = ODG;

        fprintf('  ODG:  %5.2f\n',ODGs(signum,algo,gapnum))
        fprintf('  PSM:  %5.2f\n',PSMs(signum,algo,gapnum))
        fprintf('  PSMt: %5.2f\n',PSMts(signum,algo,gapnum))
    end
end

end % gapnum

end % signum