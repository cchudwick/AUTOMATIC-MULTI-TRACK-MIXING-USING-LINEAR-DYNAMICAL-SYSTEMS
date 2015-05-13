% Harrison Zafrin and Colin Chadick
% MIR Final Project
% AUTOMATIC MULTI-TRACK MIXING USING LINEAR DYNAMICAL SYSTEMS
% -------------------------------------------------------------------------

clear;
clc;

% -------------------------------------------------------------------------
% Create fft params structure to pass to all functions 
% -------------------------------------------------------------------------

% FFT Params for Spectrogram
field1 = 'win_size'; win_size = 1028;
field2 = 'hop_size'; hop_size = win_size/2;
field3 = 'noverlap'; noverlap = win_size - hop_size;

% Create ffrparams Structure
fftparams = struct(field1, win_size,...
                   field2, hop_size,...
                   field3, noverlap);
               
% -------------------------------------------------------------------------
% Create Frame Buffering Parameters
% -------------------------------------------------------------------------

% Get sampling rate
[y, Fs] = audioread('03 California Gurls (feat. Snoop Dog.wav');

% 1 Second Frames
frame_window = Fs;

% 0.75 seconds of overlap
frame_overlap = Fs*0.75;

% -------------------------------------------------------------------------
% Weight Estimation
% -------------------------------------------------------------------------

% Import 2-Track
[x_t, fs, t] = import_audio('03 California Gurls (feat. Snoop Dog.wav');

% This length will be used to zero-pad all stems to length of track
x_len = length(x_t);

% Buffer with n_overlap
x_t_buff = buffer(x_t, frame_window, frame_overlap, 'nodelay');

% Create Stem Database
% -------------------------------------------------------------------------
filenames = {'Bass.wav'};

% Create a Cell Array filled with buffered matrices
for i=1:length(filenames)
    
    % Import Stem
    [ stem ] = import_audio( filenames{i} );

    % Make the stem and the 2-track the same length
    stem = stem(1:x_len);

    % Buffer with n_overlap
    stem_buff = buffer(stem, frame_window, frame_overlap, 'nodelay');
    
    % Put into Cell Array
    stem_data{i} = stem_buff;
end

% -------------------------------------------------------------------------
% For each frame we calculate a ground truth weight coefficient
% -------------------------------------------------------------------------
for i=1:size(x_t_buff, 2);
    
    % Grab Frame out of 2-track and compute STFT of said frame
    master_frame = x_t_buff(:,i);
    
    % Compute STFT
    master_stft = spectrogram(master_frame, fftparams.win_size, fftparams.noverlap);
    
    % Vectorize and Concatenate the STFT into one column V
    V = master_stft( : );
    
    % Get the mags    
    V = abs(V);
    
    % Now we do the same for all stems, STFT, Vectorize and Concatenate
    for j=1:length(stem_data)
        
        % Pull out stem from stem_data cell array, curr_stem is a matrix
        curr_stem = stem_data{j};
        
        % Grab the current frame in time
        stem_frame = curr_stem(:,i);
        
        % Compute RMS of the Frame to determine whether or not it's active
        rms_val = rms(stem_frame);
        
        % If the frame is active        
        if rms_val > 0.01
        
            % Compute STFT of current stem Frame
            stem_stft = spectrogram(stem_frame, fftparams.win_size, fftparams.noverlap);
            
            % Vectorize and Concatenate the STFT into one column U_k            
            U_stem = stem_stft( : );
            
            % Get the mags            
            U_stem = abs(U_stem);
            
        % Else the stem is not active in the current frame and does not contribute to V
        else
            
            % Negate the stem (THIS CANT BE RIGHT)            
            U_stem = zeros(length(V), 1);
            
        end
        
        % Load stem into combination spectra matrix
        U_matrix(:,j) = U_stem;
        
    end
    
    % Now with U matrix and V vector, NNLS to get coefficient per frame
    x = lsqnonneg(U_matrix, V);
    
    % This should contain the bass fader coefficient now as a test	 
    frame_coef(i) = x;
    
end

% Test Plot
plot(frame_coef);













