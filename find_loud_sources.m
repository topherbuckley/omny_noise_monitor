% Find the loud source, and corresponding octave band over a given
% threshold using correlation analysis
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

%% Resampling and crosscorrelations


fs_ds = 1e3; % [Hz]

% Frequency vector for the autocorrelations(pwelch gives half the spectrum)
fcorr = linspace(0,fs_ds/2-1/lH,lH)'; 



% Resample source signals and microphone signals to fs_ds = 1k
% x_ds = zeros(size(x));
% mic_ds = zeros(size(mic));

for k=1:n_sources % Loop over sources
   x_ds(:,k)=resample(x(:,k),fs_ds,fs);
end

for j=1:n_mics % Loop over microphones
   mic_ds(:,j)=resample(mic(:,j),fs_ds,fs);
end


% Design the octave band filterbank
BW = '1 octave';
N = 6;           % Filter Order
F0 = 1000;       % Center Frequency (Hz)
Fs = 48000;      % Sampling Frequency (Hz)
oneOctaveFilter = octaveFilter('FilterOrder', N, ...
    'CenterFrequency', F0, 'Bandwidth', BW, 'SampleRate', Fs)
F0 = getANSICenterFrequencies(oneOctaveFilter);
F0(F0<20) = [];
F0(F0>20e3) = [];
Nfc = length(F0);
for i=1:Nfc
    fullOctaveFilterBank{i} = octaveFilter('FilterOrder', N, ...
        'CenterFrequency', F0(i), 'Bandwidth', BW, 'SampleRate', Fs); %#ok
end






% Calculate correlations
for j=1:n_mics % Loop over microphones
    
    for k=1:n_sources % Loop over sources
        tic
%         y(k,j,:) = fftfilt(flipud(x_ds(:,k)),mic_ds(:,j));
        [y(k,j,:), tcorr] = xcorr(x_ds(:,k),mic_ds(:,j));
        Py(k,j,:) = pwelch(squeeze(y(k,j,:)),NFFT,NFFT/2,NFFT,fs_ds);
        toc
    end
end

%
figure(1); clf;
figure(2); clf;
for j=1:n_mics % Loop over microphones
    
    for k=1:n_sources % Loop over sources
       figure(1);
       subplot(n_mics,n_sources,(j-1)*n_sources+k);
       plot(tcorr,squeeze(y(k,j,:))); 
       xlim(tcorr([1 end]))
       title(['Mic',num2str(j),' Src',num2str(k)]);
       if (k==1); ylabel 'Xcorr [-]'; end
       if (j>1); xlabel 'Lags [s]'; end
       
       %
       figure(2);
       subplot(n_mics,n_sources,(j-1)*n_sources+k);
       plot(fcorr,db(squeeze(Py(k,j,:))));
       xlim([0 fs_ds/2]);
%        myspectrogram(squeeze(y(k,j,:)), fs, [18 1], @hanning, 1024, [-60 0] );
       title(['Mic',num2str(j),' Src',num2str(k)]);
       if (j>1); xlabel 'Frequency [Hz]'; end
       if (k==1); ylabel 'PSD [dB/Hz]'; end
    end
end

