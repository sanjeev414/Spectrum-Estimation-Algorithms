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


N = 512;  N1=N^2;% Number of Samples
numSignals = 3;                          % Number of incoming signals
SNR = -15;                                     % SNR values in db
omg = [0.310 0.315  0.145];          % Frequencies
M = 1024;  M1 = M+N; I =  eye(N);                                                    % Extended Number of Samples for high resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = 0:N-1;                                       % Samples

% Steering Vector Computation
eomg = linspace(0,1,M)';                 % Define Frequencies for elongated Steering Vector
eomg1=linspace(0,1,M1)';



% Signal Addition
signal = zeros(1,N);
for i = 1:numSignals
    Amp = 5*exp(1i*2*pi*rand([1,1]));
    z = Amp*exp(1i*2*pi*omg(i)*t);
    signal = signal+z;
end


% Generate noise
noise = (randn(1,N ) + 1i * randn(1,N )) * sqrt(0.5);  % Generation of Pseudo-Random Noise
noisePower = 10^(-SNR/10); % Noise Power
noise = sqrt(noisePower) * noise;
E = sqrt(noisePower).*eye(N);

% Adding Signal and Noise
signals = (signal +noise)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Power Computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



DFT = exp(-1i*2*pi*(0:N-1)'*(0:M-1)/M);
DFT = [DFT I];
IDFT = DFT';



% Find the intial power values
s_old = abs(IDFT*signals).^2/norm(DFT)^4;

%% Algorithm Implementation


for ii = 1:1e6
    R = DFT*diag(s_old)*IDFT + E;
    Rinv = inv(R);
    aa = (abs(IDFT*Rinv*signals).^2);
    bb = [];
    for ij = 1: M1
        bb = [bb;(abs(IDFT(ij,:)*Rinv*DFT(:,ij)).^2)];
    end
    s_new = aa./bb;
    Er = norm(s_new-s_old)/norm(s_old);
    if(Er>1e-4)
        s_old = s_new;
    else
        break;
    end
end

s_new=real(s_new);
snew1 = s_new(1:M);


figure; plot(eomg,snew1);


