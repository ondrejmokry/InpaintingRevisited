%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      INPAINTING METHODS COMPARISON                      %
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

global conj_atoms;
conj_atoms = true;

%% paths and settings
addpath('reweighted l1 relaxation');
addpath('SPAIN');
addpath('Janssen');
addpath('OMP');

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
turnon = logical([ 1,... %  1. DR
                   1,... %  2. CP
                   1,... %  3. wDR
                   1,... %  4. wCP
                   1,... %  5. reDR
                   1,... %  6. reCP
                   0,... %  7. gradual
                   1,... %  8. tdc
                   1,... %  9. SSPAIN H
                   0,... % 10. SSPAIN OMP
                   1,... % 11. ASPAIN
                   1,... % 12. OMP
                   1,... % 13. Janssen
                   ]);

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

prompt = '\nChoose signal number (1-10). Vector input accepted: ';
signums = input(prompt);

prompt = '\nChoose gap length number (1-10). Vector input accepted: ';
gapnums = input(prompt);

prompt = '\nChoose number of gaps per one signal (recommended max. 10): ';
n = input(prompt);

fprintf('\nControl the algorithms to be used:\n')
prompt = '  Use defaults? (0/1): ';
defaults = input(prompt);

if ~defaults
    methods = {...
        'Douglas-Rachford',...
        'Chambolle-Pock',...
        'weighted Douglas-Rachford',...
        'weighted Chambolle-Pock',...
        'iteratively reweighted Douglas-Rachford',...
        'iteratively reweighted Chambolle-Pock',...
        'gradual',...
        'direct time-domain compensation (tdc)',...
        'S-SPAIN H',...
        'S-SPAIN OMP',...
        'A-SPAIN',...
        'OMP',...
        'Janssen'};
    for method = 1:length(methods)
        prompt = sprintf('    %s (0/1): ',methods{method});
        turnon(method) = logical(input(prompt));
    end
end    

prompt = '\nPlot spectrograms? (0/1): ';
spectrograms = input(prompt);

prompt = '\nSave wavs? (0/1): ';
wavs = input(prompt);

sig_counter = 0;

%% inpainting
for signum = signums
    
sig_counter = sig_counter + 1;
gap_counter = 0;    

for gapnum = gapnums

gap_counter = gap_counter + 1;    
    
fprintf('\nSignal number: %d (no. %d of the chosen %d)',signum,sig_counter,length(signums));
fprintf('\nGap length: %d ms (no. %d of the chosen %d)\n\n',gaps(gapnum),gap_counter,length(gapnums));
    
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

% frame for methods based on convex relaxation
F = frametight(frame('dgt',{'hann',w},a,M,'timeinv'));

% parameter of the transform for methods based on non-convex heuristics
DFTred = 4;

%% signal degradation
gap_length  = gaps(gapnum);
h           = round(fs*gap_length/1000); % gap length in samples 
full.length = length(signal);
full.mask   = true(full.length,1);

%% segment setting
not_degraded   = 0.5; % length of reliable part at the start / end of the signal (is seconds)
segment.length = round((length(signal)-2*not_degraded*fs)/n);

%% some global parameters 
off           = 'half';        % 'full' / 'half' / 'none'
gradual.type  = 'analysis';    % 'analysis' / 'synthesis'
TDCparam.type = 'analysis';    % 'analysis' / 'synthesis'

%% some not so global parameters
DRweighting   = 'norm';        % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'
CPweighting   = 'energy';      % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'
TDCweighting  = 'energy';      % 'none'/ 'supp' / 'abs' / 'norm' / 'energy'

%% initializing solutions
solution.DR          = signal;
solution.CP          = signal;
solution.wDR         = signal;
solution.wCP         = signal;
solution.reDR        = signal;
solution.reCP        = signal;
solution.gradual     = signal;
solution.tdc         = signal;
solution.SSPAIN_H    = signal;
solution.SSPAIN_OMP  = signal;
solution.ASPAIN      = signal;
solution.OMP         = signal;
solution.Janssen     = signal;

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
%                            Douglas-Rachford                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(1)
        fprintf('Starting Douglas-Rachford...\n')
        
        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.DR(idxs) = segment.solution(1:segment.length)*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Chambolle-Pock                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(2)
        fprintf('Starting Chambolle-Pock...\n');
        
        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.CP(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       weighted Douglas-Rachford                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(3)
        fprintf('Starting weighted Douglas-Rachford...\n')

        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = DRweighting;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.wDR(idxs) = segment.solution(1:segment.length)*segment.max;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                        weighted Chambolle-Pock                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(4)
        fprintf('Starting weighted Chambolle-Pock...\n');
        
        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = CPweighting;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.wCP(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                iteratively reweighted Douglas-Rachford                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(5)
        fprintf('Starting iteratively reweighted Douglas-Rachford...\n');
        
        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;

        % transform settings
        param.type = 'synthesis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = true;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], RWparamsolver);

        % updating the global solution
        solution.reDR(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 iteratively reweighted Chambolle-Pock                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    if turnon(6)
        fprintf('Starting iteratively reweighted Chambolle-Pock...\n');
        
        % parameters of the main cycle
        RWparamsolver.maxit = 10;
        RWparamsolver.epsilon = 1e-3;
        RWparamsolver.delta = 0.01;
        
        % transform settings
        param.type = 'analysis';
        param.F = F;
        param.offset = off;   
        param.weighting = 'none';
        param.reweighting = true;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], RWparamsolver);

        % updating the global solution
        solution.reCP(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                gradual                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(7) && strcmp(gradual.type,'analysis')
        fprintf('Starting gradual...\n');
        
        gapindexes = find(~ssegment.mask);
        gradual.s = gapindexes(1);
        gradual.f = gapindexes(end);
        gradual.r = ceil(h/4);
        gradual.mask = ssegment.mask;
        gradual.gapped = ssegment.gapped;
        gradual.atoms = syn_atoms(F);

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
            gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, CPweighting);

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
    end  
    
    if turnon(7) && strcmp(gradual.type,'synthesis')
        fprintf('Starting gradual...\n');
        
        gapindexes = find(~ssegment.mask);
        gradual.s = gapindexes(1);
        gradual.f = gapindexes(end);
        gradual.r = ceil(h/4);
        gradual.mask = ssegment.mask;
        gradual.gapped = ssegment.gapped;
        gradual.atoms = syn_atoms(F);

        % not changing parameters        
        DRparam.g = @(x) 0; % it is actually the indicator function...
        c = frana(F,ssegment.gapped);
        DRparam.dim = length(c);
        
        % solver settings
        DRparamsolver = [];
  
        % the main cycle
        while gradual.s <= gradual.f
            % computing the weights for l1 relaxation
            gradual.wts = weights( F, L, gradual.mask, U, V, gradual.atoms, DRweighting);

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
    
    if turnon(8)      
        fprintf('Starting direct time-domain compensation (tdc)...\n');

        % (1) taking the basic reconstruction of the original gap
        
        % transform settings
        param.type = TDCparam.type;
        param.F = F;
        param.offset = off;   
        param.weighting = TDCweighting;
        param.reweighting = false;

        % algorithm
        segment.solution = reweighted(segment.gapped, segment.mask, param, [], [], []);

        % updating the global solution
        solution.tdc(idxs) = segment.solution(1:segment.length)*segment.max;
        
        % (2) direct time domain compensation (working with the whole signal)
        
        % ensuring only the single gap
        TDCmask = true(length(solution.tdc),1);
        TDCmask(idxs) = segment.mask;
        
        % TDC parameters
        TDCparam.F = F;
        TDCparam.offset = off;
        TDCparam.weighting = TDCweighting;
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

    if turnon(9)
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
%                              S-SPAIN OMP                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(10)
        fprintf('Starting S-SPAIN OMP...\n');
        
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
        SPAINparamsolver.f_update = 'OMP';

        fprintf('OMP chosen to compute the f-update. The algorithm may take up to several minutes.\n');

        % algorithm    
        [recovered_SSPAIN, ~, ~] = spain_segmentation(ssegment.gapped,SPAINparam,SPAINparamsolver,ssegment.data);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_SSPAIN(1:origL);

        % updating the global solution
        solution.SSPAIN_OMP(idxs) = segment.solution(1:segment.length)*segment.max;      
    end    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                A-SPAIN                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(11)
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
%                                  OMP                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if turnon(12)
        fprintf('Starting OMP...\n');
        
        OMPparam.w = w;
        OMPparam.a = a;
        OMPparam.wtype = 'hann';
        OMPparam.M = M;
        OMPparam.Ls = L;
        OMPparam.mask = ssegment.mask;
        OMPparam.F = frame('dft');
        OMPparam.F.redundancy = DFTred;
        OMPparam.F.frana = @(insig)dft([insig; zeros(length(insig)*(OMPparam.F.redundancy-1),1)]);
        OMPparam.F.frsyn = @(insig)postpad(idft(insig),length(insig)/OMPparam.F.redundancy);

        % solver settings
        OMPparamsolver.crit = 'mse';    % termination function
        OMPparamsolver.epsilon = 0.1/w;  % stopping criterion of termination function
        OMPparamsolver.maxit = OMPparam.w*OMPparam.F.redundancy;
        OMPparamsolver.store_snr = 0; 
        OMPparamsolver.store_obj = 0;
        OMPparamsolver.sol = 'qr'; % algorithm to compute projection in OMP
        
        % algorithm    
        recovered_OMP = omp_segmentation(ssegment.gapped,OMPparam,OMPparamsolver);

        % updating the segment solution
        segment.solution = segment.gapped;
        segment.solution(q:q+origL-1) = recovered_OMP(1:origL);

        % updating the global solution
        solution.OMP(idxs) = segment.solution(1:segment.length)*segment.max;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                Janssen                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if turnon(13)
        fprintf('Starting Janssen...\n');
        
        Janssenparam.w = w;
        Janssenparam.a = a;
        Janssenparam.wtype = 'hann';
        Janssenparam.Ls = L;
        Janssenparam.mask = ssegment.mask;

        % solver settings
        Janssenparamsolver.Nit = 50; % number of iterations
        % Janssenparamsolver.p = 2*a;    
    
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
if strcmp(gradual.type,'analysis')
    gradualweighting = CPweighting;
else
    gradualweighting = DRweighting;
end
leg = {'original',...
       'syn. $$\ell_1$$ (\texttt{none})',...
       'ana. $$\ell_1$$ (\texttt{none})',...
       ['syn. $$\ell_1$$ (\texttt{', DRweighting, '})'],...
       ['ana. $$\ell_1$$ (\texttt{', CPweighting, '})'],...
       'syn. $$\ell_1$$ (\texttt{iterative})',...   
       'ana. $$\ell_1$$ (\texttt{iterative})',... 
       [gradual.type(1:3), '. $$\ell_1$$ (\texttt{', gradualweighting '}) + gradual'],...
       [TDCparam.type(1:3), '. $$\ell_1$$ (\texttt{', TDCweighting '}) + tdc'],...
       'S-SPAIN H',...
       'S-SPAIN OMP',...
       'A-SPAIN',...
       'OMP',...
       'Janssen'};

fields = fieldnames(solution);
fields = fields(turnon);
   
figure()
hold on
plot(fulltime,signal)
for i = 1:sum(turnon)
    plot(fulltime,solution.(fields{i}));
end
legend(leg([true turnon]),'interpreter','latex')
xlabel('time [s]')
title(sprintf('signal: %s, gap length: %d ms, full signal',signame,gap_length),'interpreter','none')

% only the gaps
figure()
hold on
plot(signal(~full.mask))
for i = 1:sum(turnon)
    plot(solution.(fields{i})(~full.mask));
end
legend(leg([true turnon]),'interpreter','latex')
title(sprintf('signal: %s, gap length: %d ms, only the gaps',signame,gap_length),'interpreter','none')

%% SNRs (and ODGs)
fprintf('\nSNRs and ODGs computed from all the gaps at once:\n')
for i = 1:sum(turnon)
    restoredsignal = solution.(fields{i});
    fprintf('   %11s: SNR: %5.2f dB\n',fields{i},snr_n(signal(~full.mask),restoredsignal(~full.mask)))
    if PQ
        [~, ~, ODG, ~] = audioqual(signal, restoredsignal, fs);
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
    figure
    subplot(A,B,1)
    sg(signal)
    title(titles{1},'interpreter','latex')
    for spec = 2:number
        subplot(A,B,spec)
        sg(solution.(fields{spec-1}))
        title(titles{spec},'interpreter','latex')
    end
end

%% wavs
if wavs
    for i = 1:sum(turnon)
        restoredsignal = solution.(fields{i});
        audiowrite(['wavs\',sigs{signum}(1:3), '_', num2str(gaps(gapnum)), '_', fields{i}, '.wav'],restoredsignal,fs);
    end
    audiowrite(['wavs\',sigs{signum}(1:3), '_reference.wav'],signal,fs);
    audiowrite(['wavs\',sigs{signum}(1:3), '_', num2str(gaps(gapnum)), '_anchor.wav'],signal.*full.mask,fs);
end

end % gapnum

end % signum