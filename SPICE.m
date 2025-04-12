
%%
 % Copyright (c) 2025, Sanjeeva Reddy S
 % All rights reserved.
 
 %This source code is licensed under the MIT license found in the
 % LICENSE file in the root directory of this source tree.
 
 % Unauthorized copying of this file, via any medium, is strictly prohibited
 % unless explicit permission is granted by the copyright owner.
 
 % Description:
 % This file contains utility functions for processing sparse arrays.
 
 % Author: Sanjeeva Reddy S
 % EMail: sanjeevareddy.s414@gmail.com
 % Created on: January 5, 2025




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear the cache
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;close;tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



N = 200;  N1=N^2; N2=sqrt(N);% N is Number of Samples
numSignals = 3;                          % Number of incoming signals
SNR = 5;                                     % SNR values in db
omg = [ 0.3110 0.3150  0.1450];          % Frequencies
M = 512;   I=eye(N);  R1=M+N;                                                  % Extended Number of Samples for high resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = 0:N-1;
eomg = linspace(0,1,M)';
eomg1=linspace(0,1,R1)';
avg_s = zeros(1,M+N);
% Signal Addition

signal = zeros(1,N);
for i = 1:numSignals
    Amp = 10*exp(1i*2*pi*rand([1,1]));
    z = Amp*exp(1i*2*pi*omg(i)*t);
    signal = signal+z;
end


signals = awgn(signal,SNR,'measured');
% % Generate noise
% noise = (randn(1,N ) + 1i * randn(1,N )) * sqrt(0.5);  % Generation of Pseudo-Random Noise
% noisePower = 10^(-SNR/10); % Noise Power
% noise = sqrt(noisePower) * noise;

% Adding Signal and Noise
% signals = signal +noise;


% Finding the Normalized Power
Signal_Power = (abs(fft(signals,M))).^2;  % Computation of FFT
S_P = (Signal_Power)/N1;  % Dividing by norm of Steering Vector
S_Po=[S_P ((abs(signals).^2))];



%% SPICE Iterative Algorithm

D1 = exp(-1i*2*pi*(0:N-1)'*(0:M-1)/M);
D=[D1 I];
ID = D';


N_S=norm(signals);

w=[(N2/N_S)*ones(1,M) (1/N_S)*ones(1,(R1-M))];

s_old= S_Po';
for jj=1:1e6

    R = D*diag(s_old)*ID;
    R_inv=inv(R);
    u = R_inv*signals';
    c = abs(ID*u);
    
    rho = (w*(s_old.*c));

    s_new =(s_old.*c).*(1./(rho*w'));


    Er_En = (norm(s_old-s_new))/(norm(s_old));

    if(Er_En>1e-5)

        s_old = s_new;
    else
        break;
    end
end


s_new=real(s_new);
snew1 = s_new(1:M);


figure; plot(eomg,S_Po(1:M));
title("Power Spectrum Using Periodogram");
xlabel("Normalized Frequency");ylabel("Amplitude");

figure; plot(eomg,snew1);
title("Power Spectrum Using SPICE Optimization");
xlabel("Normalized Frequency");ylabel("Amplitude");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
