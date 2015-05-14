% Harrison Zafrin
% stem_frame = time domain frame of a signal
% win_size = size of the window
% fs = sampling rate
% -------------------------------------------------------------------------
% Compute the Spectral Centroid on a Window
% -------------------------------------------------------------------------
function [ SC ] = spectral_centroid( stem_frame, win_size, fs )

% Create fk vector
hertz_vector = linspace(0, fs, win_size);

% Get the magnitudes of the frame
stem_mags = abs(fft(stem_frame));

% Calculate the numerator
for k=1:length(stem_mags)
    numerator(k) = hertz_vector(k)*stem_mags(k);
end

% Sum across k
numerator = sum(numerator);

% Sum across k 
denominator = sum(stem_mags);

% Calc the SC for the window/frame
SC = numerator/denominator;

end

