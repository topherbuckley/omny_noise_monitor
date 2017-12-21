% Generate 1/3-octave filterbank
%
% Giuliano Bernardi
% Created:           Dec 22, 2017
% Last update:       Dec 22, 2017

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

fs = 1e3;
bw = 1/3;

fMin = 32.5;
fMax = 100;

octs = log2(fMax/fMin);
bmax = ceil(octs/bw);

fc = fMin*2.^( (0:bmax) * bw ); % centre frequencies
fl = fc*2^(-bw/2); % lower cutoffs
fu = fc*2^(+bw/2); % upper cutoffs

numBands = length(fc);

b = cell(numBands,1);
a = cell(numBands,1);

figure
for nn = 1:length(fc)

    [b{nn},a{nn}] = butter(2, [fl(nn) fu(nn)]/(fs/2), 'bandpass');
    [h,f]=freqz(b{nn},a{nn},1024,fs);

    hold on;
    plot(f, 20*log10(abs(h)) );

end
set(gca, 'XScale', 'log')
ylim([-50 0])


dt=1/fs;
t = 0:dt:2000;      % 2 seconds
fo = 0.5; f1 = 120; 
x = chirp(t,fo,numel(t)*dt,f1,'logarithmic')';
y = zeros(numel(x), length(fc));
for nn = 1:length(fc)
    y(:,nn) = filter(b{nn},a{nn},x);
end

figure; plot(y)