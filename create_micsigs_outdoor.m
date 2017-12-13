% Simulate outdoor propagation of sound. From multiple sources to multiple
% monitoring stations
%
% Giuliano Bernardi
% Created:           Dec 03, 2017
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

% Audio folder
audio_folder = 'audio';

for k = 1:11
    source_filename{k} = ['Track-',num2str(k,'%02.f'),'.wav'];
end

% Signal length
siglength = 60;


% Increase the distances increases the delay in the AIRs. 
% Therefore I have to consider what's the best choice to:
% - truncate the AIR in line 81 (so I don't loose the peak)
% - truncate the resulting microphone signal post convolution (so I don't
%   have to keep using a very long mic signal for the correlations)




% Number of monitoring mics
n_mics = 2;
% Number of stages
n_sources = 3;
% Stage levels
source_levs = [1; 0.1; 0.1];

% Sampling frequency [Hz]
fs = 16e3;

% Time vector [s]
t = (0:1/fs:siglength-1/fs);

% PSD parameters and frequency vector [Hz]
NFFT = 1024;
lH = NFFT/2+1;
f = linspace(0,fs/2-1/lH,lH)'; % Frequency vector (pwelch gives half the spectrum)



% Distance source/mic [m]
%       s1    s2    
dist = [5500 5400; ...
        5510 5405; ...
        5505 5410];

    
% Speed of sound [m/s]
c = 343; 

% Calculate delay
delay = dist/c;

% Calculate Octave Bands
fcentre = 10^3 * (2 .^ [-6:4]); % Frequencies [Hz]
fcentre(fcentre > fs/2) = [];
% Add 0 frequency (necessary for the filter design procedure)
fcentre_yw = [0 fcentre]';

% Parameters to evaluate sound propagation
T=27;    % Temperature [degree Celsius]
hr=80;   % Relative humidity in percentage hr=80 means 80 percent humidity
ps=1;    % Is the barometric pressure ratio. Usually, ps=1;

% figure(1); clf; hold on;

% Run the program
% [alpha, alpha_iso, c, c_iso]=air_absorption(f, T, hr, ps);
for j = 1:n_mics
    for k = 1:n_sources
        % Calculate IR (include inverse square law and delay)
        tmp_air = [zeros(round(delay(k,j)*fs),1); source_levs(k)/dist(k,j)^2; zeros(siglength*fs,1)];
        air(k,j,:) = tmp_air(1:siglength*fs);
        
%         % Plot the AIRs
%         plot(squeeze(air(k,j,:)));

        % Calculate attenuation values
        tmp_a = atmAtten(T,ps,hr,dist(k,j),fcentre);
        a(k,j,:) = [1; normc(1./tmp_a')];

        % Fit the filter to the octave band attenuation values
%         [byw(k,j,:),ayw(k,j,:)] = yulewalk(25,fcentre_yw/fs*2,squeeze(a(k,j,:)));
        % Using FIR2 (frequency sampling-based FIR filter design)
        bfir2(k,j,:) = fir2(1000,fcentre_yw/fs*2,squeeze(a(k,j,:)));


            
%         figure(1); clf;
%         freqz(squeeze(byw(k,j,:)),squeeze(ayw(k,j,:)),1024,fs)
%         subplot(2,1,1); hold on;
%         plot(fcentre_yw,db(squeeze(a(k,j,:))),'r')
%         
%         
%         figure(2); clf;
%         freqz(squeeze(bfir2(k,j,:)),1,1024,fs)
%         subplot(2,1,1); hold on;
%         plot(fcentre_yw,db(squeeze(a(k,j,:))),'r')
%         tilefigs([2 2])
    end
end

%%

% Other methods. Don't work very well. Didn't check much why.

% byw = firpm(100,fcentre/fs*2,a)
% byw = firls(100,fcentre/fs*2,a);


%%

% Load different speech signals
for k=1:n_sources
       [source,fs_wav]=audioread(fullfile(audio_folder,source_filename{k}));
       x(:,k)=normc(resample(source(1:siglength*fs_wav,1),fs,fs_wav));
end

   
% Create mic signals
figure(1);
for j=1:n_mics % Loop over microphones
   mic(:,j)=zeros(size(x,1),1);
   
   subplot(n_mics,1,j); cla; hold on;
   
   for k=1:n_sources % Loop over sources
       tic
%        src_mic_contribution = filter(squeeze(byw(k,j,:)),squeeze(ayw(k,j,:)),fftfilt(squeeze(air(k,j,:)),x(:,k)));
%        src_mic_contribution = fftfilt(squeeze(bfir2(k,j,:)),fftfilt(squeeze(air(k,j,:)),x(:,k)));
       src_mic_contribution = fftfilt(squeeze(air(k,j,:)),x(:,k));
       plot(t,src_mic_contribution);
       mic(:,j) = mic(:,j) + src_mic_contribution;
       % Here I should add environmental noise
       toc
   end
   title(num2str(dist(:,j)'));
end

save('mic','mic','fs')