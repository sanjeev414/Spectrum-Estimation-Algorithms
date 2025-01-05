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
SNR = -15;                                     % SNR values in db
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


signals = awgn(signal,SNR);
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
    % w=[];rho=0;
    % for k2=1:R1
    %     ww=norm(D(:,k2))./norm(signals);
    %     w=[w ww];
    %     rho=rho+(ww*s_old(k2)*c(k2));
    % end
    rho = (w*(s_old.*c));

    s_new = (1/rho)*(s_old.*c).*(1./w');


    Er_En = (norm(s_old-s_new))/(norm(s_old));

    if(Er_En>1e-5)

        s_old = s_new;
    else
        break;
    end
end
s_new=real(s_new);
snew1 = s_new(1:M);

% figure;plot(eomg1,(s_new));
figure; plot(eomg,S_P);
title("Power Spectrum Using Periodogram");
xlabel("Normalized Frequency");ylabel("Amplitude");

figure; plot(eomg,snew1);
title("Power Spectrum Using SPICE Optimization");
xlabel("Normalized Frequency");ylabel("Amplitude");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% SNR estimation
% 
% p = sort(s_new,'ascend');
% sn = zeros(M,1);
% for n0 = 0:M-1
%     for k0 = 0:M-1
%         sn(n0+1) = p(k0+1)/ (n0+1);
%     end
% end
% 
% 
% qn = zeros(M,1);
% for n1 = 0:M-1
%     for k1 = 0:M-1
%         qn(n1+1) = (((p(k1+1)^2)/(n1+1))-sn(n1+1)^2);
%     end
% end
% 
% tn = (sn.^2)./qn;
% 
% k1 = [];k2 = [];
% for ii=0:M-1
%     if(tn(ii+1)<1)
%         k1 = [k1,(ii+1)];
%     elseif (tn(ii+1)>=1)
%         k2 = [k2,(ii+1)];
%     end
% end
% 
% ps = sum(p(k1))/length(k1);
% 
% pn = sum(p(k2))/length(k2);
% 
% est_SNR = ps/pn;




