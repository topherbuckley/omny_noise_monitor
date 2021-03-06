% Find the loud source, and corresponding octave band over a given
% threshold using correlation analysis
%
% Giuliano Bernardi
% Created:           Dec 13, 2017
% Last update:       Dec 25, 2017

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

%% Resampling and crosscorrelations

clear all; close all; clc;

% Add path of support functions
addpath(genpath('support_functions'));

load(fullfile(tempdir,'computed_scenario')); % Load computed scenario

fs_ds = 1e3; % Resampled fs [Hz]

% Threshold
th = -240;


% PSD parameters and frequency vector [Hz]
NFFT = 1024;
lH = NFFT/2+1;
f = linspace(0,fs/2-1/lH,lH)'; % Frequency vector (pwelch gives half the spectrum)

% Frequency vector for the autocorrelations(pwelch gives half the spectrum)
fcorr = linspace(0,fs_ds/2-1/lH,lH)'; 


% Create third-octave filter bank
T = 100e-3; % Integration time, i.e. window length
[B,A] = adsgn(fs_ds); 



% Resample source signals and microphone signals to fs_ds = 1k
% x_ds = zeros(size(x));
% mic_ds = zeros(size(mic));

for k=1:n_sources % Loop over sources
   x_ds(:,k)=resample(x(:,k),fs_ds,fs);
end

for j=1:n_mics % Loop over microphones
   mic_ds(:,j)=resample(mic(:,j),fs_ds,fs);
end


% Initialize max values
max_y = 0;
max_Py = 0;
max_Poct3 = -inf;

% Calculate correlations
for j=1:n_mics % Loop over microphones
    
    for k=1:n_sources % Loop over sources
        tic
%         y(k,j,:) = fftfilt(flipud(x_ds(:,k)),mic_ds(:,j));
        [y(k,j,:), tcorr(k,j,:)] = xcorr(x_ds(:,k),mic_ds(:,j));
        [~, imax(k,j)] = max(abs(y(k,j,:)));
        sel_inds(k,j,:) = imax(k,j)-5*fs_ds:imax(k,j)+5*fs_ds-1;
        y_tr(k,j,:) = y(k,j,squeeze(sel_inds(k,j,:)));
                
        Py(k,j,:) = pwelch(squeeze(y_tr(k,j,:)),NFFT,NFFT/2,NFFT,fs_ds);
        [tmp_Poct3, tmp_Foct3] = filtbank(squeeze(y_tr(k,j,:)),fs_ds,T,'extended');
        Poct3{k,j} = tmp_Poct3; Foct3{k,j} = tmp_Foct3;
        exc_inds{k,j} = find(mean(Poct3{k,j})>th == 1);
        toc
        
        % Find max values
        max_y = max([max_y abs(max(y_tr(k,j,:)))]);
        max_Py = max([max_Py max(Py(k,j,:))]);
        max_Poct3 = max([max_Poct3 max(mean(Poct3{k,j}))]);
        
    end
end

%%
figure(1); clf;
figure(2); clf;
figure(3); clf;
figure(4); clf;
for j=1:n_mics % Loop over microphones
    
    for k=1:n_sources % Loop over sources
       curr_tcorr = tcorr(squeeze(sel_inds(k,j,:)))/fs_ds; 
       curr_tcorr_Poct3 = curr_tcorr(1:T*fs_ds:end);
        
       figure(1);
       subplot(n_mics,n_sources,(j-1)*n_sources+k);
       plot(curr_tcorr,squeeze(y_tr(k,j,:))); 
       ylim(max_y*[-1.15 1.15])
       title(['Mic',num2str(j),' Src',num2str(k)]);
       if (k==1); ylabel 'Xcorr [-]'; end
       if (j>1); xlabel 'Lags [s]'; end
       
       %
       figure(2);
       subplot(n_mics,n_sources,(j-1)*n_sources+k);
       plot(fcorr,db(squeeze(Py(k,j,:))));
       xlim([fcorr(1) fs_ds/2]);
       ylim([db(max_Py)-150 db(max_Py)]);
       set(gca,'xscale','log');
       set(gca,'XTick',[1 10 100 1000]);
       title(['Mic',num2str(j),' Src',num2str(k)]);
       if (j>1); xlabel 'Frequency [Hz]'; end
       if (k==1); ylabel 'PSD [dB/Hz]'; end
       
       %
       figure(3);
       subplot(n_mics,n_sources,(j-1)*n_sources+k);
       imagesc(Poct3{k,j}'); axis xy
       set(gca,'YTick',[2:3:length(Foct3{k,j})]);
       set(gca,'YTickLabel',Foct3{k,j}(2:3:length(Foct3{k,j})));
       set(gca,'XTick',[20:20:length(curr_tcorr_Poct3)]);
       set(gca,'XTickLabel',round(10*curr_tcorr_Poct3(1:length(curr_tcorr_Poct3)/5:end))/10);
       title(['Mic',num2str(j),' Src',num2str(k)]);
       if (j>1); xlabel 'Time [s]'; end
       if (k==1); ylabel 'Frequency [Hz]'; end       
       
       %
       figure(4);
       mean_Poct3 = mean(Poct3{k,j},1);
       subplot(n_mics,n_sources,(j-1)*n_sources+k); hold on;
       plot(Foct3{k,j},mean_Poct3); axis xy
       title(['Mic',num2str(j),' Src',num2str(k)]);
       xlim([Foct3{k,j}(1) Foct3{k,j}(end)]);
       ylim([(max_Poct3)-50 (max_Poct3)]);
       set(gca,'XTick',round(Foct3{k,j}(1:2:end)));
       set(gca,'XTickLabel',round(Foct3{k,j}(1:2:end)));
       set(gca,'xscale','log'); 
       for i = 1:length(exc_inds{k,j})
           stem(Foct3{k,j}(exc_inds{k,j}),mean_Poct3(exc_inds{k,j}),'o')
       end
       if (j>1); xlabel 'Frequency [Hz]'; end
       if (k==1); ylabel 'PSD [dB/oct]'; end       
       
    end
end

