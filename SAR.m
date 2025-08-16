    %% SAR Simulation with Cube Target
    
    % Physical constants
    c = physconst('LightSpeed'); % Speed of light (m/s)
    fc = 10e9;                  % Radar center frequency (Hz, X-band)
    lambda = c / fc;            % Wavelength (m)
    
    % System parameters
    v = 150;                    % Platform velocity (m/s)
    h = 1000;                   % Platform altitude (m)
    prf = 1000;                 % Pulse Repetition Frequency (Hz)
    tpd = 3e-6;                 % Pulse width (s)
    bw = 150e6;                 % Signal bandwidth (Hz)
    fs = 2 * bw;                % Sampling frequency (Hz)
    
    % Radar platform definition
    radarPlatform = phased.Platform('InitialPosition', [0; 0; h],'Velocity', [0; v; 0]);
    
    %% Cube target parameters
    cube_center = [1000; 0; 0];         % position (m)
    cube_dims = [10; 10; 10];           % size [Lx; Ly; Lz] (m)
    points_per_side = 6;                % scatter points per cube edge
    rcs_per_point = 0.5;                 % RCS of each scatter point (m²)
    
    % Generate cube scatter points
    [targetPositions, targetRCS] = generateCubeTargets(cube_center, cube_dims, points_per_side, rcs_per_point);
    
    % Define radar target and platform
    cubeTarget = phased.RadarTarget('OperatingFrequency', fc,'MeanRCS', targetRCS);
    cubePlatform = phased.Platform('InitialPosition', targetPositions,'Velocity', zeros(3, size(targetPositions, 2)));
    
    %% Waveform (Linear FM)
    waveform = phased.LinearFMWaveform('SampleRate', fs, 'PulseWidth', tpd,'PRF', prf, 'SweepBandwidth', bw, 'OutputFormat', 'Pulses');
    
    % Antenna
    antenna = phased.CosineAntennaElement;
    antennaGain = aperture2gain(1, lambda); % Isotropic with 1 m² aperture
    radiator = phased.Radiator('Sensor', antenna, 'OperatingFrequency', fc, 'PropagationSpeed', c);
    collector = phased.Collector('Sensor', antenna, 'PropagationSpeed', c, 'OperatingFrequency', fc);
    
    % Transmitter and Receiver
    transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
    receiver = phased.ReceiverPreamp('SampleRate', fs, 'NoiseFigure', 30);
    
    % Propagation Channel
    channel = phased.FreeSpace('PropagationSpeed', c, 'OperatingFrequency', fc, 'SampleRate', fs, 'TwoWayPropagation', true);
    
    %% Simulation parameters
    apertureDuration = 4; % Total time for synthetic aperture (s)
    numpulses = floor(apertureDuration * prf);
    
    % Number of fast-time samples per pulse
    num_fasttime_samples = round(waveform.SampleRate * waveform.PulseWidth);
    
    % Preallocate raw data (fast-time × slow-time)
    raw_data = zeros(num_fasttime_samples, numpulses);
    
    %% Main simulation loop
    t1 = tic;
    
    for m = 1:numpulses
        [radarPos, radarVel] = radarPlatform(1/prf);
        [cubePos, cubeVel] = cubePlatform(1/prf);
    
        rx_total = zeros(num_fasttime_samples, 1);
        t2 = tic;
        for k = 1:size(cubePos, 2)
            [range, angle] = rangeangle(cubePos(:,k), radarPos);
    
            tx_signal = transmitter(waveform());
            tx_signal = radiator(tx_signal, angle);
            rx_signal = channel(tx_signal, radarPos, cubePos(:,k), radarVel, cubeVel(:,k));
    
            % Create a target with the RCS of just this point
            ptTarget = phased.RadarTarget('OperatingFrequency', fc, 'MeanRCS', targetRCS(k));
    
            rx_signal = ptTarget(rx_signal);
            rx_signal = collector(rx_signal, angle);
            rx_signal = receiver(rx_signal);
    
            rx_total = rx_total + rx_signal(1:num_fasttime_samples);
        end
        
        elapsed = toc(t2);  % stores the elapsed time in seconds
        fprintf('Inside loop %d took %.3f seconds.\n', k, elapsed);
        raw_data(:, m) = rx_total;
        elapsed = toc(t1);  % stores the elapsed time in seconds
        fprintf('Outside Loop %d took %.3f seconds.\n',m,elapsed);
    end
    
    %% Range Compression
    matchedFilter = phased.MatchedFilter('Coefficients', getMatchedFilter(waveform));
    range_compressed = matchedFilter(raw_data);
    
    %% Azimuth Compression
    [num_range_bins, num_pulses] = size(range_compressed);
    pulse_index = 0:(num_pulses-1);
    fd = (2 * v / lambda) * sin(0); % Broadside Doppler
    
    azimuth_ref = exp(-1j * 2 * pi * fd * pulse_index / prf);
    azimuth_compressed = range_compressed .* repmat(azimuth_ref, num_range_bins, 1);
    
    %% SAR Image Formation
    sar_image = fftshift(fft(azimuth_compressed, [], 2), 2);
    
    %% Display
    figure;
    imagesc(abs(sar_image));
    colormap(gray);
    title('Final SAR Image of Cube Target');
    xlabel('Azimuth Bins');
    ylabel('Range Bins');
    

