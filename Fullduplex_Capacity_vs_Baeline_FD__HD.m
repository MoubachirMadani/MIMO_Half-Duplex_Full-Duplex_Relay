%%% Please cite this article "Half-Duplex and Full-Duplex Interference
%%% Mitigation in Relays Assisted Heterogeneous Network" PlosOne ID: PONE-D-22-34295
%by Moubachir Madani Fadoul and Chee-Onn Chow"
%%% Figure number 3.
clearvars
close all
clc

% number of channel realization
It = 1000;
Ps=1; Pr=1;
% Number of transmit and receive antennas (2-4-6-6)or (2-6-2-4) or
% (2-2-4-6)- antt(2-4-6-4)
M=2;  % no. of antenas at source / destination 
Nr=4;  % receive antenna at relay
Mr=6;  % transmit antenna at the relay
Nd=4;  % transmit antenna at the relay
% SNR range in dB
% SNR = ratio between total transmit power and noise variance
SNRdBvalues =  [0:2:20]; %[-10:2:30];
Czf     = zeros(1,It);
Cmmse   = zeros(1,It);

SNRidx = 0;
for SNRdB=SNRdBvalues
    SNRdB
    SNRidx = SNRidx + 1;
    SNR=10^(SNRdB/10);
    
    % compute maximal achievable rate 
    for kk=1:It
        
        % generate channel realization
        HS1R1 = (randn(Nr,M)+j*randn(Nr,M))/sqrt(2);
        HR1D1 = ( randn(Nd,Mr) + j*randn(Nd,Mr))/sqrt(2);
        HR1 = ( randn(Nr,Mr) + j*randn(Nr,Mr))/sqrt(2);
        HR1R2 = ( randn(Nr,Mr) + j*randn(Nr,Mr))/sqrt(2);
        HR1D2 = ( randn(Nd,Mr) + j*randn(Nd,Mr))/sqrt(2);
         
        HS2R2 = ( randn(Nr,M) + j*randn(Nr,M))/sqrt(2);        
        HR2D2 = ( randn(Nd,Mr) + j*randn(Nd,Mr))/sqrt(2);
        HR2 = ( randn(Nr,Mr) + j*randn(Nr,Mr))/sqrt(2);
        HR2R1 = ( randn(Nr,Mr) + j*randn(Nr,Mr))/sqrt(2);
        HR2D1 = ( randn(Nd,Mr) + j*randn(Nd,Mr))/sqrt(2);  
         [U S V] = svd(HS1R1);
        [U1 S1 V1] = svd(HR1D1);
 %% Null at the relay tx 
  Wtu2=null(HR1R2); 
  Wt2=Wtu2';
  
 PSHRR2 =pinv(HR1);  
Wtu1=(PSHRR2*HR1R2*Wt2');  
Wt1=Wtu1';

% Null at the destination Rx 
 Wrd1=null(HR2D1); 
 Wrd2=null(HR1D2); 
  
%% Prpoposed Scheme -(S1-R1)i
            UU=(SNR/M *HS1R1*HS1R1');
%               RR=SNR/Mr *(Wr1'*HR1'*HR1*Wr1);  %RSI ->becomes zero
%               YY=SNR/Mr *(Wr1'*HR2D1'*Wt2'*Wt2*HR2D1*Wr1');  %RDI  ->zeros          
           %  CS1R1(kk)=real(log2(det(eye(M)+(AA/RR+YY))));  %(1:2,1:2)
         CS1R1(kk)=real(log2(det(eye(Nr)+(UU))));  
%  Prpoposed Scheme of-(Ri-Di) 
         PP=SNR/Mr* Wrd1'*HR1D1'*HR1D1*Wrd1;
%          YY=SNR/Mr *(Wr1'*HR2D1*Wt2'*Wt2*HR2D1'*Wr1);  %RDI  ->zeros 
        CR1D1(kk)=real(log2(det(eye(M)+(PP))));
     
        %% Prpoposed Scheme -(S2-R2)j
            AA=(SNR/M *HS2R2*HS2R2');
            DR=SNR/Mr *(HR1*Wt2'*Wt2*HR1');  %RSI ->becomes zero
%               YY=SNR/Mr *(Wr1'*HR2D1'*Wt2'*Wt2*HR2D1*Wr1');  %RDI  ->zeros          
           %  CS1R1(kk)=real(log2(det(eye(M)+(AA/RR+YY))));  %(1:2,1:2)
         CS2R2(kk)=real(log2(det(eye(Nr)+(AA))));  
         
%  Prpoposed Scheme of-(R-D) 
         DP=SNR/Mr* Wrd2'*HR2D1'*HR2D1*Wrd2;
%          YY=SNR/Mr *(Wr1'*HR2D1*Wt2'*Wt2*HR2D1'*Wr1);  %RDI  ->zeros 
        CR2D2(kk)=real(log2(det(eye(M)+(DP))));
          end
            
     % min
        CAPS1RD1=min(CS1R1,CR1D1); %Proposed Scheme -eq.(15)
        CAPS2RD2=min(CS2R2,CR2D2); %Proposed Scheme -eq.(15)
        CAPSRD=min(CAPS1RD1,CAPS2RD2); 

    % average over all channel realizations for a given SNR value
    CAPS1R1D1(SNRidx)  = mean(CAPSRD);
        
end

% plot
figure(1)

 plot(SNRdBvalues, CAPS1R1D1,'k-+','linewidth',2.3, 'MarkerSize',10)
 
title('Capacity of SR-hop and RD-hops')
xlabel(' SNR (dB)') %\rho
ylabel('Capacity (bits/sec/Hz)')
 legend('Proposed Scheme-equal PA')
