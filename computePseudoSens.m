function [absSinf, Swz, Swy] = computePseudoSens(freqs, sysR, A_rho, frfGwz, frfGuz, frfGwy, frfGuy, nrHOSIDFsMaxUser, SamplesHighestHarmonicUser)
    % ---------------------------------------------------------------------
    % This function can be utilized to compute the pseudo-sensitivity 
    % magnitude (from w to z) of a Lure-type reset control system:
    %
    %             z   _____    w 
    %           <----|     |<-----
    %             y  |  G  |   u
    %           -----|_____|<-----
    %           |     _____      |
    %           |    |     |     |
    %           |--->|  R  |-----|
    %                |_____|
    % 
    % Signal z and the underlying higher-order sinusoidal-input sensitivity
    % functions (HOSISFs) are computed using [1, Theorem 3.1]. The pseudo-
    % sensitivity is computed based on the maximum of z, as in (35) of [1].
    % 
    % v00 - Luke van Eijk (15/03/2025)
    % Code based on:
    % [1] L.F. van Eijk, D. Kostić, S.H. HosseinNia, "Frequency Response Analysis
    %       of Lure-Type Reset Control Systems," submitted to IEEE Control Systems Letters
    % ---------------------------------------------------------------------
    % Input definition:
    %   freqs - 1-by-M linearly-spaced frequency array (Hz), such that
    %               freqs = [f_1, f_2, ..., f_M] with f_k = k*f_1
    %   sysR  - struct containing variables A_R, B_R, and C_R of system R
    %               as in (1) of [1]
    %   A_rho - variable A_rho of system R as in (1) of [1]
    %   frfGwz - 1-by-M complex-valued array with frequency-response function
    %               (FRF) of SISO LTI element Gwz at frequencies 'freqs'
    %   frfGuz - 1-by-M complex-valued array with FRF of Guz at frequencies 'freqs'
    %   frfGwy - 1-by-M complex-valued array with FRF of Gwy at frequencies 'freqs'
    %   frfGuy - 1-by-M complex-valued array with FRF of Guy at frequencies 'freqs'
    %   nrHOSIDFsMaxUser (optional) - Largest HOSIDF that should be taken into account
    %   SamplesHighestHarmonicUser (optional) - Number of to be evaluated time-instants
    %               per time-period of highest-harmonic (see [1, Section V])
    %
    % Output definition:
    %   absSinf - 1-by-M array with pseudo-sensitivity magnitudes at frequencies 'freqs'
    %   Swz     - M-by-M (or nrHOSIDFsMax-by-M) array with higher-order sinusoidal-input
    %               sensitivity functions (HOSISFs) at frequencies 'freqs' (from w to z)
    %   Swy     - M-by-M (or nrHOSIDFsMax-by-M) array with HOSISFs at frequencies 'freqs' (from w to y)
    % ---------------------------------------------------------------------

    nrFreqs = length(freqs);    % Number of frequencies
    if nargin == 7
        nrHOSIDFs = nrFreqs;            % Maximum number of HOSIDFs that can be taken into account
        SamplesHighestHarmonic = 100;   % By default, a very accurate but computationally expensive setting
    elseif nargin == 8
        nrHOSIDFs = min(nrHOSIDFsMaxUser,nrFreqs);  % Maximum can (optionally) be lowered using input 'nrHOSIDFsMaxUser'
        SamplesHighestHarmonic = 100;
    elseif nargin == 9
        nrHOSIDFs = min(nrHOSIDFsMaxUser,nrFreqs);
        SamplesHighestHarmonic = SamplesHighestHarmonicUser;  % Accuracy/computation-time trade-off can (optionally) be changed
    else
        error('An unexpected error occured.')
    end


    %% Compute HOSISFs    
    % Compute HOSIDFs of the reset element in (1) of [1]
    A_R = sysR.A_R; B_R = sysR.B_R; C_R = sysR.C_R;
    [Hn, Swy, Swz] = deal(NaN(nrHOSIDFs,nrFreqs));
    for nn = 1:nrHOSIDFs
        Hn(nn,:) = computeResetHOSIDF(A_R, B_R, C_R, 0, A_rho, freqs, nn); % compute HOSIDFs
    end
    frfRbl = computeResetHOSIDF(A_R, B_R, C_R, 0, eye(length(B_R)), freqs, 1); % FRF of reset-element's base-linear system: (2) in [1]

    % Compute HOSISFs of reset control system in Fig. 4 of [1]
    for nn = 1:nrHOSIDFs
        if nn == 1
            Swy(1,:) = frfGwy ./ (1 - frfGuy .* Hn(1,:));      % (16) in [1]
            Swz(1,:) = frfGwz + frfGuz .* Hn(1,:) .* Swy(1,:); % (19) in [1]
        else
            omegaIdxs = 1:floor(nrFreqs/nn);
            nnOmegaIdxs = nn * omegaIdxs;

            frfDummy = Hn(nn,omegaIdxs) .* Swy(1,omegaIdxs) .* exp(1i*(nn-1)*angle(Swy(1,omegaIdxs))) ./ ...
                                                                                        (1 - frfGuy(nnOmegaIdxs) .* frfRbl(nnOmegaIdxs));
            Swy(nn,omegaIdxs) = frfGuy(nnOmegaIdxs) .* frfDummy;    % (17) in [1]
            Swz(nn,omegaIdxs) = frfGuz(nnOmegaIdxs) .* frfDummy;    % (20) in [1]
        end

        % figure(101)
        % semilogx(freqs,mag2db(abs(Swz(nn,:))))
        % hold on
        % 
        % figure(102)
        % semilogx(freqs,mag2db(abs(Swy(nn,:))))
        % hold on
    end  


    %% Compute pseudo-sensitivity magnitudes
    absSinf = NaN(1,nrFreqs);
    for kk = 1:nrFreqs
        n_max = min(nrHOSIDFs,floor(nrFreqs/kk)); % Largest HOSIDF that can be taken into account for this input frequency, (33) in [1]

        freqInput = freqs(kk);          % Input frequency (Hz)
        omegaInput = 2*pi*freqInput;    % Input frequency (rad/s)
        Tperiod = 1 / freqInput;        % Period of external input (s)
        nrSamples = SamplesHighestHarmonic * n_max;     % Number of time-instants to evaluate
        time = Tperiod * ((1:nrSamples)-1) / nrSamples; % (34) in [1]
        perfOutput = zeros(1,nrSamples); %resetInput = zeros(1,nrSamples);
        for nn = 1:n_max
            % Without loss of generality, we assume \hat{w} = 1 & \varphi_w = 0 (see [1, Section V])
            perfOutput = perfOutput + abs(Swz(nn,kk)) * sin(nn*omegaInput*time + angle(Swz(nn,kk)));    % Performance output based on (18) in [1]
            % resetInput = resetInput + abs(Swy(nn,kk)) * sin(nn*omegaInput*time + angle(Swy(nn,kk)));    % Reset element's input based on (15) in [1]
        end
        absSinf(kk) = max(abs(perfOutput));    % Pseudo-sensitivity magnitude as in (35) of [1]

        % figure
        % plot(time,perfOutput)    
    end
end