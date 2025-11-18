% lab20_ex_IQpoints.m
% Text conversion to carrier states numbers and (Ik,Qk) components, and back
clear all; close all;

text_in = 'Hello World! 01234567890 abcdefgh',
modtype = '64QAM';     % 2PAM, 4PAM, 8PAM, BPSK, QPSK, DQPSK, 8PSK, 4QAM, 16QAM
do_texttransmit = 0;  % 0/1 optional further processing
do_decode = 0;        % 0/1 optional text decoding

% Definition of constellation points of carrier states
[IQcodes, Nstates, Nbits, R ] = IQdef( modtype );
    phi = 0:pi/100:2*pi; c=R*cos(phi); s=R*sin(phi);
    figure; plot(c,s,'k-',real(IQcodes), imag(IQcodes),'ro','MarkerFaceColor','red'); 
    xlabel('I(k)'); ylabel('Q(k)'); title('Possible IQ(k) states'); grid; pause

% Coding our message using carrier states
numbers = text2numbers( text_in, Nbits ), pause    % text to IQ state numbers
%numbers  = floor( Nstates*(rand(100,1)-10*eps) ); % random states
IQk = numbers2IQ( numbers, modtype, IQcodes );     % IQ state numbers to IQ values
    figure;
    subplot(211); stem(real(IQk),'b'); grid; xlabel('k'); title('I(k)');
    subplot(212); stem(imag(IQk),'r'); grid; xlabel('k'); title('Q(k)'); pause

% Data transmission - our next exercise
% IMPORTANT: set "if(0)" in the beginning of the program lab20_ex_pulse_shaping.m
if( do_texttransmit ) lab20_ex_pulse_shaping, end   % our next exercise

% Text decoding - our next exercise - carrier state decoding
if( do_decode == 1 )
    numbers = IQ2numbers( IQk, modtype )             % find carrier state numbers 
    text_out = numbers2text( numbers, Nbits ), pause % convert them to text
end


% B = log2(N) - N = constelation points