function [data_rec_fin]=omp_segmentation(data_gapped, param, paramsolver)
% OMP_SEGMENTATION serves as a middle layer between the main files and the
% OMP algorithm
% This fucntion performs input signal padding, computation of analysis and
% synthesis window, windowing the signal and after processing each block by
% OMP it folds the blocks back together.
%
% Input parameters
%       data_gapped    vector of gapped signal
%       param          structure of parameters containing:
%                            Ls         length of the original signal
%                            w          window length (in samples)
%                            a          window shift (in samples)
%                            wtype      string defining the type of the window, see http://ltfat.github.io/doc/sigproc/firwin.html
%                            F          frame (usually DFT frame)
%                            mask       logical vector indicating the reliable samples
%
%       paramsolver    structure of parameters for OMP containing:
%                            sol        solver of the projection during OMP
%                            epsilon    stopping threshold of termination function
%                            maxit      maximal possible number of iterations with particular settings
%                            store_snr  switch to enable computing SNR in each iteration
%                            store_obj  switch to enable storing the value of termination function in each iteration
%                            conj_atoms switch to enable choosing of pairs of complex conjugate atoms in each iteration of OMP
%
% Output parameters
%       data_rec_fin   vector of reconstructed signal after OLA
%
% Date: 08/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% preparation for padding at the end of the signal
L = ceil(param.Ls/param.a)*param.a+(ceil(param.w/param.a)-1)*param.a; % L is divisible by a and minimum amount of zeros equals gl (window length). Zeros will be appended to avoid periodization of nonzero samples.
N = L/param.a; % number of signal blocks

% padding the signals and mask to length L
padding = zeros(L-param.Ls, 1);
data_gapped = [data_gapped; padding];
param.mask = [param.mask; true(L-param.Ls,1)];

% construction of analysis and synthesis windows
g = gabwin(param.wtype, param.a, param.w, L);
gana = normalize(g,'peak'); % peak-normalization of the analysis window
gsyn = gabdual(gana, param.a, param.w)*param.w; % computing the synthesis window

% this is substituting fftshift (computing indexes to swap left and right half of the windows)
idxrange = [0:ceil(param.w/2)-1,-floor(param.w/2):-1];
idxrange2 = idxrange+abs(min(idxrange))+1;

% initialization of param_seg (parameters for one signal block)
param_seg = param;
param_seg.Ls = param.w;
param_seg.mask = true(param.w,1);

% initialization of signal blocks
data_block = zeros(param.w,1);
data_rec_fin = zeros(L,1);

% dictionary of DFT
full_dict = dft(eye(param.w*param.F.redundancy))';
full_dict = full_dict(1:param.w,:);

for n=0:N-1       
    % multiplying signal block with windows and choosing corresponding masks
    idx = mod(n*param.a + idxrange,L) + 1;
    data_block(idxrange2) = data_gapped(idx).*gana;
    param_seg.mask(idxrange2) = param.mask(idx);
    
    % if the block does not contain any samples to restore, skip it
    if sum(param_seg.mask) == param.w
        continue;
    end
    
    % dictionary reduction
    dict = full_dict(param_seg.mask,:);
    [~, width] = size(dict);
    
    % vector for reduced dictionary normalization
    w = zeros(width,1);
    for i = 1:width
        w(i) = 1/norm(dict(:,i));
    end
    
    % construction of P and Pt for greed_omp_qr(...)
    % for functions P, Pt, see end of this file

    % OMP  
    [x, ~, ~] = greed_omp_qr(data_block(param_seg.mask),@P,width,...
        'P_trans',@Pt,...
        'stopCrit',paramsolver.crit,...
        'stopTol',paramsolver.epsilon*sum(param_seg.mask)/param.w,...
        'maxIter',paramsolver.maxit);
    data_rec_block = frsyn(param.F,w.*x);
    data_rec_block = real(data_rec_block);
    
    % Folding blocks together using Overlap-add approach (OLA)
    data_rec_block = ifftshift(data_rec_block);
    data_rec_fin(idx) = data_rec_fin(idx) + data_rec_block.*gsyn;
end

% ensure equality of solution and data in reliable part
data_rec_fin(param.mask) = data_gapped(param.mask);

% crop the padding of reconstructed signal
data_rec_fin = data_rec_fin(1:param.Ls);

function syn = P(x)
    syn = frsyn(param.F,w.*x);
    syn = syn(param_seg.mask);
end

function ana = Pt(y)
    z = zeros(param.w,1);
    z(param_seg.mask) = y;
    ana = w.*frana(param.F,z);
end

end

