function [frfGwz, frfGuz, frfGwy, frfGuy] = convertToLure(frfC1, frfC2, frfC3, frfC4, frfC5, frfPlant)
    % ---------------------------------------------------------------------
    % This function converts the reset control system below into the Lure
    % form (as in Fig. 4 of [1]) with reference 'r' as external input 'w'
    % and error 'e' as performance output 'z'
    %
    %          --> C2 --> R --> C3 --
    %          |                    |                                      
    % --> C1 --|                    + --> C5 --> P -->
    %          |                    |                
    %          ---------> C4 --------
    % 
    %                   |
    %                   |
    %                  \|/
    %            z    _____    w 
    %           <----|     |<-----
    %            y   |  G  |   u
    %           -----|_____|<-----
    %           |     _____      |
    %           |    |     |     |
    %           |--->|  R  |-----|
    %                |_____|
    %
    % v00 - Luke van Eijk (15/03/2025)
    % Code based on:
    % [1] L.F. van Eijk, D. Kostić, S.H. HosseinNia, "Frequency Response Analysis
    %       of Lure-Type Reset Control Systems," submitted to IEEE Control Systems Letters
    % ---------------------------------------------------------------------
    % Input definition:
    %   frfC1 - 1-by-M (or M-by-1) complex-valued array with frequency-response function (FRF)
    %                                                            of SISO LTI controller C1
    %   frfC2 - 1-by-M (or M-by-1) complex-valued array with FRF of SISO LTI controller C2
    %   frfC3 - 1-by-M (or M-by-1) complex-valued array with FRF of SISO LTI controller C3
    %   frfC4 - 1-by-M (or M-by-1) complex-valued array with FRF of SISO LTI controller C4
    %   frfC5 - 1-by-M (or M-by-1) complex-valued array with FRF of SISO LTI controller C5
    %   frfPlant - 1-by-M (or M-by-1) complex-valued array with FRF of SISO LTI plant P
    %
    % Output definition:
    %   frfGwz - 1-by-M (or M-by-1) complex-valued array with FRF of SISO
    %               LTI element from external input 'w' to performance output 'z'
    %   frfGuz - 1-by-M (or M-by-1) complex-valued array with FRF of SISO
    %               LTI element from 'u' to 'z'
    %   frfGwy - 1-by-M (or M-by-1) complex-valued array with FRF of SISO
    %               LTI element from 'w' to 'y'
    %   frfGuy - 1-by-M (or M-by-1) complex-valued array with FRF of SISO
    %               LTI element from 'u' to 'y'
    % ---------------------------------------------------------------------

    %% Convert to form in Fig. 3 of [1]
    %            ---> R -----
    %            |          |                
    % --> Cpre --|          + --> Cpos --> P -->
    %            |          |                
    %            --> Cpar ---

    frfCpre = frfC1 .* frfC2;
    frfCpar = frfC4 ./ frfC2 ./ frfC3;
    frfCpos = frfC3 .* frfC5;


    %% Convert to Lure-form in Fig. 4 of [1]
    frfGwz = 1 ./ (1 + frfPlant .* frfCpos .* frfCpar .* frfCpre);  % (21) in [1]
    frfGuz = -frfPlant .* frfCpos .* frfGwz;                        % (22) in [1]
    frfGwy = frfCpre .* frfGwz;                                     % (23) in [1]
    frfGuy = frfCpre .* frfGuz;                                     % (24) in [1]
end