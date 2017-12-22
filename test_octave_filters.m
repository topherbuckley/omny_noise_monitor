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

clear all; close all; clc;

addpath(genpath('support_functions'));


% Time specifications:
fs_ds = 1e3;                 % [Hz]
dt = 1/fs_ds;                % [s/sample]
StopTime = 5;                % [s]
t_tr = (0:dt:StopTime-dt)';  % [s]
fcs = [63 120 250];          % [Hz]
x_tr = zeros(size(t_tr));   
% Sine waves
for k=1:length(fcs)
    x_tr = x_tr + cos(2*pi*fcs(k)*t_tr);
end

% Frequency specifications
% PSD parameters and frequency vector [Hz]
NFFT = 1024;
lH = NFFT/2+1;
f = linspace(0,fs_ds/2-1/lH,lH)'; % Frequency vector (pwelch gives half the spectrum)


% Create third-octave filter bank
T = 100e-3; % Integration time, i.e. window length
t_P = (0:T:StopTime-T)';  % [s]
[B,A] = adsgn(fs_ds); 
% x_trfb = filter(B,A,x_tr); 
[P,F] = filtbank(x_tr,fs_ds,T,'extended');


%%
% Spectrum with welch's method
Px = pwelch(x_tr,NFFT,NFFT/2,NFFT,fs_ds);

% Spectrum with third-octave filters
T_long = StopTime; % Integration time equal to StopTime
[P_long,F_long] = filtbank(x_tr,fs_ds,T_long,'extended');  % Using long T
[P_oct3,F_oct3] = oct3bank(resample(x_tr,44.1e3,fs_ds));   % Using oct3bank 


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
logfsgram(x_tr,NFFT,fs_ds,NFFT,NFFT/2,15);
title 'Log-frequency spectrogram'
subplot(4,1,4);
% bankdisp(P,F,-40,-20);
imagesc(t_P,[1:length(F)],P'); axis xy
set(gca,'YTick',[2:3:length(F)]);
set(gca,'YTickLabel',F(2:3:length(F)));
title 'Octave-band spectrogram'
%%
figure(2); clf;
subplot(2,1,1);
plot(f,db(Px));
xlim([0 max(F_long)]);
subplot(2,1,2); hold on;
plot(F,mean(P,1));
plot(F,mean(P_long,1));
% plot(F,mean(P_short,1));
plot(F_oct3,P_oct3);
xlim([0 max(F_long)]);

% imagesc(P'); axis xy
% waterfall(P);
% zlabel('Level [dB]');
% ylabel('Frame #');
% xlabel('Frequency band');
% set(gca,'XTick',[2:3:length(F)]);
% set(gca,'XTickLabel',F(2:3:length(F)));