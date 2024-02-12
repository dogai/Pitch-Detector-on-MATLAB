% CMPE468 - Term Project: Pitch Detector

while true
    disp('Enter the audio file name (including the file extension), or enter "exit" to quit:');
    file_name = input('', 's');
    file_name = strtrim(file_name); % Remove leading/trailing whitespace
    
    if strcmp(file_name, 'exit')
        disp('Exiting the program...');
        break;
    end
    
    try

        % Step 1: Read the audio file
        [y, fs] = audioread(file_name);

        y_normalized = y / max(abs(y)); % Normalize the signal by dividing it by the maximum absolute value

        cutoff_frequency = 4000;
        normalized_cutoff = cutoff_frequency / (fs / 2);
        [b, a] = butter(6, normalized_cutoff, 'low');
        y_filtered = filter(b, a, y_normalized);

        frame_duration = 30e-3;
        frame_length = round(frame_duration * fs);
        num_frames = floor(length(y_filtered) / frame_length);
        frames = reshape(y_filtered(1:num_frames * frame_length), frame_length, num_frames);
        spectra = abs(fft(frames));
        spectral_mean = mean(spectra, 2);
        noise_threshold = spectral_mean * 1.5;
        spectra_denoised = max(spectra - noise_threshold, 0);
        frames_denoised = ifft(spectra_denoised, 'symmetric');
        y_denoised = reshape(frames_denoised, [], 1);


        % Step 2: Frame Segmentation
        frame_duration = 30e-3;
        frame_length = round(frame_duration * fs);
        num_frames = floor(length(y_denoised) / frame_length);
        frames = reshape(y_denoised(1:num_frames * frame_length), frame_length, num_frames);


        % Step 3: Spectral Segmentation
        window = hamming(frame_length);
        frames_windowed = frames .* window;
        spectra = abs(fft(frames_windowed));


        % Step 4: Pitch Estimation (Autocorrelation method)
        pitch_estimates = zeros(1, num_frames);

        for i = 1:num_frames
            frame = frames(:, i);
            autocorr = xcorr(frame);
            [~, max_lag] = max(autocorr);

            pitch_estimates(i) = fs / (max_lag - 1);
        end


        % Step 5: Pitch Refinement
        refined_pitch_estimates = pitch_estimates;

        for i = 1:num_frames
            frame = frames(:, i);
            spectrum = abs(fft(frame));
            hps = spectrum;
            for k = 2:4
                downsampled_frame = downsample(frame, k);
                upsampled_spectrum = interpft(abs(fft(downsampled_frame)), length(frame));
                hps = hps .* upsampled_spectrum;
            end
            [~, max_index] = max(hps);
            refined_pitch_estimates(i) = fs / max_index;
        end


        % Step 6: Pitch Tracking and Analysis
        pitch_track = zeros(1, length(y));

        for i = 1:num_frames
            frame_start = (i - 1) * frame_length + 1;
            frame_end = frame_start + frame_length - 1;
            pitch_track(frame_start:frame_end) = refined_pitch_estimates(i);
        end


% Step 7: Visualization and Output
        time = (0:length(y)-1) / fs; % Create a time axis for visualization

        % Plot the waveform and overlay the pitch track
        figure;
        subplot(2,1,1);
        plot(time, y);
        xlabel('Time (s)');
        ylabel('Amplitude');
        title('Waveform');

        subplot(2,1,2);
        plot(time, pitch_track);
        xlabel('Time (s)');
        ylabel('Pitch (Hz)');
        title('Pitch Track');

        % Save the pitch track as a text file
        output_file = 'pitch_track.txt';
        dlmwrite(output_file, pitch_track, 'delimiter', '\t');

        % Display a message indicating that the pitch track is saved
        disp('Pitch track saved as "pitch_track.txt"');
        drawnow; % Flush the output buffer and update the display
    catch
        disp('Error: Unable to read the audio file. Please make sure the file exists and is in the correct format.');
        drawnow; % Flush the output buffer and update the display
    end
end