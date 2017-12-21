% Test octave filters and compare to simple spectrogram and log-spaced
% spectrogram
%
% Giuliano Bernardi
% Created:           Dec 13, 2017
% Last update:       Dec 13, 2017

%
% This code is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your %
% option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see http://www.gnu.org/licenses/.

addpath(genpath('support_functions'));


% Time specifications:
fs_ds = 1e3;                 % [Hz]
dt = 1/fs_ds;                % [s/sample]
StopTime = 1;                % [s]
t_tr = (0:dt:StopTime-dt)';  % [s]
fcs = [63 120 250];          % [Hz]
x_tr = zeros(size(t_tr));   
% Sine waves
for k=1:length(fcs)
    x_tr = x_tr + cos(2*pi*fcs(k)*t_tr);
end


% Create third-octave filter bank
T = 100e-3; % Integration time, i.e. window length
[B,A] = adsgn(fs_ds); 
% x_trfb = filter(B,A,x_tr); 
[P,F] = filtbank(x_tr,fs_ds,T,'extended');


% [S_sp,F_sp,T_sp] = spectrogram(x_tr,512,1024,512,fs,'yaxis');


figure(1); clf;
subplot(4,1,1);
plot(t_tr,x_tr);
title 'Time signal'
subplot(4,1,2);
myspectrogram(x_tr, fs_ds, [18 1], @hanning, 1024, [-60 0] );
ylim([0 fs_ds/2]);
title 'Spectrogram'
subplot(4,1,3);
logfsgram(x_tr,1024,fs_ds);
title 'Log-frequency spectrogram'
subplot(4,1,4);
% bankdisp(P,F,-40,-20);
imagesc(P'); axis xy
set(gca,'YTick',[2:3:length(F)]);
set(gca,'YTickLabel',F(2:3:length(F)));
title 'Octave-band spectrogram'




% imagesc(P'); axis xy
% waterfall(P);
% zlabel('Level [dB]');
% ylabel('Frame #');
% xlabel('Frequency band');
% set(gca,'XTick',[2:3:length(F)]);
% set(gca,'XTickLabel',F(2:3:length(F)));