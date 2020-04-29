function y = inpaintFrame_janssenInterpolation(problemData,param)
% Frame-level inpainting method based on the linear prediction by
% Janssen.
%
% Usage: xEst = inpaintFrame_janssenInterpolation(problemData,param)
%
%
% Inputs:
%          - problemData.x - observed signal to be inpainted
%          - problemData.Imiss - Indices of clean samples
%          - param.p - Order of the autoregressive model used in
%          for linear prediction
%
% Outputs:
%          - y: estimated frame
%
% References:
% 
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).

s = problemData.x;
N = length(s);
Im = find(problemData.IMiss);
IObs = find(~problemData.IMiss);
M = length(Im);
Im = sort(Im); Im = Im(:); % Im: indexes of missing samples
s(Im) = 0;

if nargin<2 || ~isfield(param,'p')
   p = min(3*M+2,round(N/3));
else
   p = param.p;
end
if nargin<2 || ~isfield(param,'GR')
   param.GR = false;
end
if nargin<2 || ~isfield(param,'NIt')
   NIt = 100;
else
   NIt = param.NIt;
end


%IAA = abs(Im*ones(1,N)-ones(M,1)*(1:N));
IAA = abs(repmat(Im,1,N)-repmat(1:N,M,1));
IAA1 = IAA<=p;
AA = zeros(size(IAA));

if param.GR
   figure;
   hold on
end
for k=1:NIt
   
   % Re-estimation of LPC
   aEst = lpc(s,p).';
   
   % Re-estimation of the missing samples
   b = aEst.'*hankel(aEst.',[aEst(end),zeros(1,p)]);
   AA(IAA1) = b(IAA(IAA1)+1);
   % xEst = -inv(AA(:,Im))*AA(:,IObs)*s(IObs); % use Chol to invert matrix
   [R, flagErr] = chol(AA(:,Im));
   if flagErr
      % xEst = - AA(:,Im)\(AA(:,IObs)*s(IObs));
      xEst = -inv(AA(:,Im))*AA(:,IObs)*s(IObs);
   else
      xEst = -R\(R'\(AA(:,IObs)*s(IObs)));
   end
   s(Im) = xEst;
   if param.GR
      e = filter(aEst,1,s);
      plot(k,10*log10(mean(e(p+1:end).^2)),'o')
   end
end

y = s;

return
