% Simulate outdoor propagation of sound. From multiple sources to multiple
% monitoring stations
%
% Giuliano Bernardi
% Created:           Dec 03, 2014
% Last update:       Dec 06, 2017

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

source_filename{1}='speech1.wav';
source_filename{2}='speech2.wav';
source_filename{3}='Babble_noise1.wav';

siglength=20;


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
source_levs = [1; 0.2; 0.1];

% Sampling frequency [Hz]
fs = 16e3;

% Distance source/mic [m]
%       s1    s2    
dist = [5500 200; ...
        1010 180; ...
        2000 190];

    
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

% Run the program
% [alpha, alpha_iso, c, c_iso]=air_absorption(f, T, hr, ps);
for j = 1:n_mics
    for k = 1:n_sources
        % Calculate IR (include inverse square law and delay)
        tmp_air = [zeros(round(delay(k)*fs),1); source_levs(k)/dist(k,j)^2; zeros(siglength*fs,1)];
        air(k,j,:) = tmp_air(1:siglength*fs);

        % Calculate attenuation values
        tmp_a = atmAtten(T,ps,hr,dist(k,j),fcentre);
        a(k,j,:) = [1; normc(1./tmp_a')];

        % -----------------------
        % -----------------------
        % -----------------------
        % I could use the function fir2 in this case instead of yulewalk. 
        % I think it's easier to use
        % -----------------------
        % -----------------------
        % -----------------------

        % Fit the filter to the octave band attenuation values
        [bf(k,j,:),af(k,j,:)] = yulewalk(25,fcentre_yw/fs*2,squeeze(a(k,j,:)));


%         figure(1); clf;
%         freqz(squeeze(bf(k,j,:)),squeeze(af(k,j,:)),1024,fs)
%         subplot(2,1,1); hold on;
%         plot(fcentre_yw,db(squeeze(a(k,j,:))),'r')
    end
end

%%

% Other methods. Don't work very well. Didn't check much why.

% bf = firpm(100,fcentre/fs*2,a)
% bf = firls(100,fcentre/fs*2,a);


%%

% Load different speech signals
clear x
for k=1:n_sources
       [source,fs_wav]=audioread(source_filename{k});
       x(:,k)=normc(resample(source(1:siglength*fs_wav,1),fs,fs_wav));
end
   
% Create mic signals
figure(1);
for j=1:n_mics % Loop over microphones
   mic(:,j)=zeros(size(x,1),1);
   
   subplot(n_mics,1,j); cla; hold on;
   
   for k=1:n_sources % Loop over sources
       tic
       src_mic_contribution = filter(squeeze(bf(k,j,:)),squeeze(af(k,j,:)),fftfilt(squeeze(air(k,j,:)),x(:,k)));
       plot(src_mic_contribution);
       mic(:,j) = mic(:,j) + src_mic_contribution;
       % Here I should add environmental noise
       toc
   end
   title(num2str(dist(:,j)'));
end

save('mic','mic','fs')

% -----------------------------------------


% corr_length=4*fs_RIR;
return
%%
for j=1:n_mics % Loop over microphones
    
    for k=1:n_sources % Loop over sources
        y(k,j,:) = fftfilt(flipud(x(:,k)),mic(:,j));
    end
end

%%
figure(1); clf;
for j=1:n_mics % Loop over microphones
    
    for k=1:n_sources % Loop over sources
       subplot(n_mics,n_sources,(j-1)*n_sources+k);
       plot(squeeze(y(k,j,:))); title(['Mic',num2str(j),' Src',num2str(k)]);
    end
end

