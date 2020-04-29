function  [snr_val] = snr_n (u, v)
% SNR_N computes computes SNR (signal-to-noise ratio) between two signals
% for evaluating quality of the inpainting restoration
%
% Input parameters
%       u     vector representing the original clean signal
%       v     vector representing the processed (reconstructed) signal  

snr_val = 20*log10(norm(u)/norm(u-v));

end