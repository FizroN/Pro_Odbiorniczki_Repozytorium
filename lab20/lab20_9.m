% lab20_ex_IQpoints.m
% This script is for Exercise 20.9, combining Listings 20.1, 20.2, and 20.5.
% It tests the conversion of text -> numbers -> IQ symbols -> numbers -> text.
% 64-QAM has been added as requested by the exercise.

clear all; close all;

% --- Parameters for Exercise 20.9 ---
text_in = 'Hello World! 0123456789 abcdefghijkl';

% --- Change this line to test different modulations ---
modtype = '64QAM'; % 2PAM, 4PAM, 8PAM, BPSK, QPSK, 8PSK, 4QAM, 16QAM, 64QAM

% As per Ex 20.9: Set do_texttransmit=0 (skip pulse shaping)
do_texttransmit = 0; 

% As per Ex 20.9: Set do_decode=1 (enable text decoding)
do_decode = 1;

% --- Main script from Listing 20.2 ---
fprintf('Testing modulation: %s\n', modtype);

% 1. Definition of constellation points
[IQcodes, Nstates, Nbits, R] = IQdef(modtype);

figure;
plot(real(IQcodes), imag(IQcodes), 'ro', 'MarkerFaceColor', 'red');
title(['Possible IQ(k) states for ', modtype]);
xlabel('I(k)'); ylabel('Q(k)'); grid on; axis equal;
if(Nstates > 16)
    % Add labels for complex QAM grids
    for k=1:Nstates
        text(real(IQcodes(k))+0.05, imag(IQcodes(k)), dec2bin(k-1, Nbits));
    end
end
pause(0.1);

% 2. Coding our message using carrier states
fprintf('Input text: %s\n', text_in);
[numbers, bitsnum, bitschar] = text2numbers(text_in, Nbits);
fprintf('Generated %d state numbers.\n', length(numbers));

% 3. Text to IQ state numbers conversion
IQk = numbers2IQ(numbers, modtype, IQcodes);

figure;
subplot(211); stem(real(IQk), 'b.'); grid; xlabel('k'); title('I(k)');
subplot(212); stem(imag(IQk), 'r.'); grid; xlabel('k'); title('Q(k)');
pause(0.1);

% 4. Data transmission (skipped)
if(do_texttransmit) 
    % lab20_ex_pulse_shaping % This part is skipped
    disp('Skipping pulse shaping (do_texttransmit = 0).');
else
    % To test decoding, we pass the "perfect" IQk signal
    % (In a real system, IQk would be corrupted by noise here)
    IQk_received = IQk; 
end

% 5. Text decoding (as per Ex 20.9)
if(do_decode == 1)
    fprintf('Decoding symbols back to text...\n');
    
    % Find carrier state numbers
    % This is the function from Listing 20.5
    numbers_out = IQ2numbers(IQk_received, modtype); 
    
    % Convert them to text
    % This is the function from Listing 20.1
    text_out = numbers2text(numbers_out, Nbits); 
    
    fprintf('Output text: %s\n', text_out);
    
    % Check for errors
    if(isequal(text_in, text_out(1:length(text_in))))
        fprintf('Success! Decoded text matches input text.\n');
    else
        fprintf('Error! Decoded text does not match.\n');
        fprintf('Input:  %s\n', text_in);
        fprintf('Output: %s\n', text_out);
    end
end


% =====================================================================
% HELPER FUNCTION: Listing 20.1 (Text <-> Numbers)
% =====================================================================

function [numbers, bitsnum, bitschar] = text2numbers(text, Nbits)
    % text to IQ state numbers conversion
    bitschar = dec2bin(double(text), 8)'; % text-array, letters in rows
    [rows, cols] = size(bitschar);
    N = rows * cols; % number of all bits
    work = reshape(bitschar, [1, N]); % bits in one row
    Nadd = Nbits - rem(N, Nbits);
    if Nadd == Nbits, Nadd = 0; end % Handle case where N is multiple of Nbits
    for k = 1:Nadd, work = [work, '0']; end % appending 0 bits at the end
    bitsnum = reshape(work, [Nbits, (N+Nadd)/Nbits])'; % bits of all states
    numbers = bin2dec(bitsnum); % state numbers: from 0 to 2^Nbits-1
end

function text = numbers2text(numbers, Nbits)
    % IQ state numbers to text conversion
    text_bits = dec2bin(numbers, Nbits);
    [rows, cols] = size(text_bits);
    text_stream = reshape(text_bits', [1, rows*cols]); % one big stream of chars '0' '1'
    N = length(text_stream); 
    N = N - rem(N, 8); % Find last full byte
    text_stream = text_stream(1:N); % remove appended bits
    text_bytes = reshape(text_stream, [8, N/8])'; % strings of bytes
    text = char(bin2dec(text_bytes))'; % conversion to text
end


% =====================================================================
% HELPER FUNCTION: Listing 20.2 (IQ Definition)
% =====================================================================

function [IQcodes, Nstates, Nbits, R] = IQdef(modtype)
    % Defines the constellation points (IQ codes) for modulations
    Nstates = 0; Nbits = 0; R = 1; IQcodes = [];
    if(isequal(modtype, '2PAM'))
        Nstates = 2; Nbits = 1; IQcodes = [-1, 1];
    elseif(isequal(modtype, '4PAM'))
        Nstates = 4; Nbits = 2; IQcodes = [-3, -1, 3, 1]; % Gray code
    elseif(isequal(modtype, '8PAM'))
        Nstates = 8; Nbits = 3; IQcodes = [-7, -5, -3, -1, 7, 5, 3, 1]; % Gray
    elseif(isequal(modtype, 'BPSK'))
        Nstates = 2; Nbits = 1; IQcodes = [1, -1]; % 0, pi
    elseif(isequal(modtype, 'QPSK'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2); 
        IQcodes = [R+j*R, -R+j*R, -R-j*R, R-j*R]; % R*exp(j*pi/4*[1 3 5 7])
    elseif(isequal(modtype, 'DQPSK'))
        Nstates = 4; Nbits = 2; R = 1;
        % Phase rotations, not absolute states
        IQcodes = [exp(j*pi/4), exp(j*3*pi/4), exp(j*7*pi/4), exp(j*5*pi/4)]; % 00, 01, 10, 11
    elseif(isequal(modtype, '8PSK'))
        Nstates = 8; Nbits = 3; 
        % Gray coded 8-PSK
        idx = [0, 1, 3, 2, 7, 6, 4, 5];
        IQcodes(idx+1) = exp(j*pi/4*(0:7));
    elseif(isequal(modtype, '4QAM'))
        Nstates = 4; Nbits = 2; R = 1/sqrt(2);
        % 00, 01, 11, 10 (Gray coded)
        IQcodes = [R*(-1-j*1), R*(-1+j*1), R*(1+j*1), R*(1-j*1)];
    elseif(isequal(modtype, '16QAM'))
        Nstates = 16; Nbits = 4; R = 1/sqrt(10);
        [I,Q] = meshgrid([-3, -1, 1, 3], [-3, -1, 1, 3]);
        IQ_map = R*(I + 1j*Q);
        % Gray code mapping
        gray_idx = [0, 1, 3, 2, 4, 5, 7, 6, 12, 13, 15, 14, 8, 9, 11, 10];
        IQcodes(gray_idx+1) = IQ_map(:);
    elseif(isequal(modtype, '64QAM'))
        % --- Added for Exercise 20.9 ---
        Nstates = 64; Nbits = 6; R = 1/sqrt(42);
        [I,Q] = meshgrid([-7,-5,-3,-1,1,3,5,7], [-7,-5,-3,-1,1,3,5,7]);
        IQ_map = R*(I + 1j*Q);
        % Standard (non-Gray) mapping for simplicity, can be improved
        IQcodes = IQ_map(:);
    else
        error('Unknown modulation type: %s', modtype);
    end
end

function IQk = numbers2IQ(numbers, modtype, IQstates)
    % State numbers to IQ values
    if(isequal(modtype, 'DQPSK'))
        IQk(1) = 1; % Initial state
        for k = 1:length(numbers)
            IQk(k+1) = IQk(k) * IQstates(numbers(k)+1); % Differential coding
        end
    else
        IQk = IQstates(numbers + 1);
    end
    IQk = IQk(:); % Ensure column vector
end


% =====================================================================
% HELPER FUNCTION: Listing 20.5 (IQ -> Numbers)
% =====================================================================

function numbers = IQ2numbers(IQ, modtype)
    % from [I,Q] values to carrier state numbers for many input IQ pairs
    N = length(IQ);
    numbers = zeros(N, 1);
    
    if(isequal(modtype, '2PAM')) % BPSK is similar
        for ns = 1:N
            if(real(IQ(ns)) < 0)
                numbers(ns) = 0; else numbers(ns) = 1;
            end
        end
    elseif(isequal(modtype, '4PAM'))
        for ns = 1:N
            I = real(IQ(ns));
            if(I < -2), numbers(ns) = 0; % -3
            elseif(I < 0), numbers(ns) = 1; % -1
            elseif(I < 2), numbers(ns) = 3; %  1
            else, numbers(ns) = 2; %  3
            end
        end
    elseif(isequal(modtype, '8PAM'))
        for ns = 1:N
            I = real(IQ(ns));
            if(I < -6), numbers(ns) = 0; % -7
            elseif(I < -4), numbers(ns) = 1; % -5
            elseif(I < -2), numbers(ns) = 3; % -3
            elseif(I < 0), numbers(ns) = 2; % -1
            elseif(I < 2), numbers(ns) = 7; %  1
            elseif(I < 4), numbers(ns) = 6; %  3
            elseif(I < 6), numbers(ns) = 5; %  5
            else, numbers(ns) = 4; %  7
            end
        end
    elseif(isequal(modtype, 'BPSK'))
         for ns = 1:N
            if(real(IQ(ns)) > 0)
                numbers(ns) = 0; else numbers(ns) = 1;
            end
         end
    elseif(isequal(modtype, '4QAM') || isequal(modtype, 'QPSK'))
        for ns = 1:N
            I = real(IQ(ns));
            Q = imag(IQ(ns));
            if(I < 0 && Q < 0), numbers(ns) = 0; % 00
            elseif(I < 0 && Q > 0), numbers(ns) = 1; % 01
            elseif(I > 0 && Q > 0), numbers(ns) = 2; % 11
            else, numbers(ns) = 3; % 10
            end
        end
    elseif(isequal(modtype, 'DQPSK'))
        % Note: This decodes the state, not the differential bits
        % A full DQPSK decoder would compare phase(k) vs phase(k-1)
        IQdiff = IQ(2:end) .* conj(IQ(1:end-1)); % Get phase change
        numbers = zeros(N-1, 1);
        for ns = 1:N-1
            ph = angle(IQdiff(ns));
            if(ph > 0 && ph < pi/2), numbers(ns)=0; % pi/4
            elseif(ph > pi/2 && ph < pi), numbers(ns)=1; % 3pi/4
            elseif(ph < -pi/2 && ph > -pi), numbers(ns)=3; % -3pi/4 (or 5pi/4)
            else, numbers(ns) = 2; % -pi/4 (or 7pi/4)
            end
        end
    elseif(isequal(modtype, '8PSK'))
        for ns = 1:N
            ph = angle(IQ(ns)); % Get phase
            ph = mod(ph, 2*pi); % Map from -pi:pi to 0:2pi
            % Find closest angle
            region = floor(ph / (pi/4)); % Quantize to 8 regions
            switch(region)
                case 0, numbers(ns) = 0; % 0
                case 1, numbers(ns) = 1; % pi/4
                case 2, numbers(ns) = 3; % pi/2
                case 3, numbers(ns) = 2; % 3pi/4
                case 4, numbers(ns) = 7; % pi
                case 5, numbers(ns) = 6; % 5pi/4
                case 6, numbers(ns) = 4; % 3pi/2
                case 7, numbers(ns) = 5; % 7pi/4
            end
        end
    elseif(isequal(modtype, '16QAM'))
        R = 1/sqrt(10);
        levels = R * [-2, 0, 2]; % Decision boundaries
        for ns = 1:N
            I = real(IQ(ns));
            Q = imag(IQ(ns));
            
            if(I < levels(1)), i_idx = 0; % -3
            elseif(I < levels(2)), i_idx = 1; % -1
            elseif(I < levels(3)), i_idx = 2; % 1
            else, i_idx = 3; % 3
            end
            
            if(Q < levels(1)), q_idx = 0; % -3
            elseif(Q < levels(2)), q_idx = 1; % -1
            elseif(Q < levels(3)), q_idx = 2; % 1
            else, q_idx = 3; % 3
            end
            
            % Map (i_idx, q_idx) back to Gray code number
            gray_map = [ 0  1  3  2;
                         4  5  7  6;
                        12 13 15 14;
                         8  9 11 10];
            numbers(ns) = gray_map(q_idx*4 + i_idx + 1);
        end
    elseif(isequal(modtype, '64QAM'))
        % --- Added for Exercise 20.9 ---
        R = 1/sqrt(42);
        % Decision boundaries are between levels
        levels = R * [-6, -4, -2, 0, 2, 4, 6];
        for ns = 1:N
            I = real(IQ(ns));
            Q = imag(IQ(ns));
            
            % Find I index (0-7)
            i_idx = 0;
            for k = 1:length(levels)
                if I > levels(k), i_idx = k; end
            end
            
            % Find Q index (0-7)
            q_idx = 0;
            for k = 1:length(levels)
                if Q > levels(k), q_idx = k; end
            end
            
            % Map (i_idx, q_idx) to number (0-63)
            % This is a simple (non-Gray) mapping
            numbers(ns) = i_idx * 8 + q_idx;
        end
    else
        disp('Demodulation is not supported!');
    end
end