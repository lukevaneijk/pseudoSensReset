function Hn = computeResetHOSIDF(A_R, B_R, C_R, D_R, A_rho, freqs, n)
    % ---------------------------------------------------------------------
    % This function computes a higher-order sinusoidal input describing
    % function (HOSIDF) of a reset controller as in (1) of [1]
    % 
    % v00 - Luke van Eijk (15/03/2025)
    % Code based on:
    % [1] N. Saikumar, K. Heinen, S.H. HosseinNia, "Loop-shaping for reset control systems:
    %       a higher-order sinusoidal-input describing functions approach," Control Engineering Practice, 2021
    % ---------------------------------------------------------------------
    % Input definition:
    %   A_R, B_R, C_R, D_R, A_rho - variables of reset controller as in (1) of [1]
    %   freqs - 1-by-M array with frequencies (Hz)
    %   n     - HOSIDF order (-)
    %
    % Output definition:
    %   Hn    - 1-by-M complex-valued array with the reset controller's 
    %               n^th-order HOSIDF, as in [1, Theorem 3.1]
    % ---------------------------------------------------------------------

    if min(freqs) <= 0
        error('Only positive frequencies are allowed')
    end
    
    % Formulas below taken from (2) in [1]
    nrFreqs = length(freqs);
    if mod(n,2) == 0 % Even orders
        Hn = zeros(1,nrFreqs);  % (13) in [1]
    elseif mod(n,2) == 1 % Odd orders
        Hn = zeros(1,nrFreqs);

        omega = 2*pi*freqs;
        nrStates = length(B_R);   
        In = eye(nrStates);
        for ff = 1:nrFreqs
            Lambda = omega(ff)^2 * In + A_R^2;                                  % (14) in [1]
            Delta = In + expm(pi*A_R/omega(ff));                                % (14) in [1]
            Delta_r = In + A_rho * expm(pi*A_R/omega(ff));                      % (14) in [1]
            Gamma_r = inv(Delta_r) * A_rho * Delta * inv(Lambda);               % (14) in [1]
            Theta_D = -2*omega(ff)^2 / pi * Delta * (Gamma_r - inv(Lambda));    % (14) in [1]

            if n == 1
                Hn(ff) = C_R * inv(1i*omega(ff)*In - A_R) * (In + 1i*Theta_D) * B_R + D_R;  % (13) in [1]
            else
                Hn(ff) = C_R * inv(1i*n*omega(ff)*In - A_R) * 1i*Theta_D * B_R; % (13) in [1] (see also corrigendum)
            end
        end
    else
        error('Only natural numbers are allowed for the HOSIDF order')
    end
end