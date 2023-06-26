%%% Please cite this article "Madani Fadoul M, Chow C-O (2023) Half-duplex and full-duplex interference mitigation in relays assisted heterogeneous network. 
%%%PLoS ONE 18(6): e0286970. https://doi.org/10.1371/journal.pone.0286970"
% Figure number 4.
clearvars
close all
clc;

% number of channel realization
It = 1000;
% Number of transmit and receive antennas
M=2;  Nr=4;   Mr=6;  Nd=4;    %(2;6;2;4)-(3;6;2;4) -(3;6;2;4);
M1=3; Nr1=6;  Mr1=9; Nd1=6;
M2=4; Nr2=8;  Mr2=12; Nd2=8;

SNRdBvalues = [0:5:20];

SNRidx = 0;
for SNRdB=SNRdBvalues
    SNRdB;
    SNRidx = SNRidx + 1;
    SNR=10^(SNRdB/10);
     for kk=1:It
        
   % generate channel realization M=2; Nr=4; Mr=6; Nd=4;
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
          [U1 S1 V1] = svd(HS1R1);
        [U11 S11 V11] = svd(HR1D1);
%% Null at the relay tx 
  %% Null at the relay tx 
  Wtu1=null(HR1); 
  Wto1=Wtu1';
  
  Wrr2=HR1R2*Wto1';
  Wrel = Wrr2';
  Wro2=null(Wrel); 
  
 Wteu1=null(HR1); 
 Wte1=Wteu1';
 
 Wteu2=null(HR2R1); 
 Wte2=Wteu2';
% Null at the destination Rx 
Wdess=HR2D1*Wte2';
 Wred1=null(Wdess'); 
 
Wdees=HR1D2*Wte1';
 Wred2=null(Wdees'); 
               
%          M1=3; Nr1=6; Mr1=9; Nd=6;
    % generate channel realization
        HS1R12 = (randn(Nr1,M1)+j*randn(Nr1,M1))/sqrt(2);
        HR1D12 = ( randn(Nd1,Mr1) + j*randn(Nd1,Mr1))/sqrt(2);
        HR12 = ( randn(Nr1,Mr1) + j*randn(Nr1,Mr1))/sqrt(2);
        HR1R22 = ( randn(Nr1,Mr1) + j*randn(Nr1,Mr1))/sqrt(2);
        HR1D22 = ( randn(Nd1,Mr1) + j*randn(Nd1,Mr1))/sqrt(2);
         
        HS2R22 = ( randn(Nr1,M1) + j*randn(Nr1,M1))/sqrt(2);        
        HR2D22 = ( randn(Nd1,Mr1) + j*randn(Nd1,Mr1))/sqrt(2);
        HR22 = ( randn(Nr1,Mr1) + j*randn(Nr1,Mr1))/sqrt(2);
        HR2R12 = ( randn(Nr1,Mr1) + j*randn(Nr1,Mr1))/sqrt(2);
        HR2D12 = ( randn(Nd1,Mr1) + j*randn(Nd1,Mr1))/sqrt(2);   
     [U2 S2 V2] = svd(HS1R12);
     [U22 S22 V22] = svd(HR1D12);
%% Null at the relay tx 
   %% Null at the relay tx 
  Wtu12=null(HR12); 
  Wto12=Wtu12';
  
  Wrr22=HR1R22*Wto12';
  Wrel2 = Wrr22';
  Wro22=null(Wrel2); 
  
 Wteu12=null(HR12); 
 Wte12=Wteu12';
 
 Wteu22=null(HR2R12); 
 Wte22=Wteu22';
% Null at the destination Rx 
Wdess2=HR2D12*Wte22';
 Wred12=null(Wdess2'); 
 
Wdees2=HR1D22*Wte12';
 Wred22=null(Wdees2'); 
        
%         M2=4; Nr2=8; Mr2=12; Nd=8;
    % generate channel realization
        HS1R13 = (randn(Nr2,M2)+j*randn(Nr2,M2))/sqrt(2);
        HR1D13 = ( randn(Nd2,Mr2) + j*randn(Nd2,Mr2))/sqrt(2);
        HR13 = ( randn(Nr2,Mr2) + j*randn(Nr2,Mr2))/sqrt(2);
        HR1R23 = ( randn(Nr2,Mr2) + j*randn(Nr2,Mr2))/sqrt(2);
        HR1D23 = ( randn(Nd2,Mr2) + j*randn(Nd2,Mr2))/sqrt(2);
         
        HS2R23 = ( randn(Nr2,M2) + j*randn(Nr2,M2))/sqrt(2);        
        HR2D23 = ( randn(Nd2,Mr2) + j*randn(Nd2,Mr2))/sqrt(2);
        HR23 = ( randn(Nr2,Mr2) + j*randn(Nr2,Mr2))/sqrt(2);
        HR2R13 = ( randn(Nr2,Mr2) + j*randn(Nr2,Mr2))/sqrt(2);
        HR2D13 = ( randn(Nd2,Mr2) + j*randn(Nd2,Mr2))/sqrt(2); 
          [U3 S3 V3] = svd(HS1R13);
        [U33 S33 V33] = svd(HR1D13);
   %% Null at the relay tx 
  Wtu13=null(HR13); 
  Wto13=Wtu13';
  
  Wrr23=HR1R23*Wto13';
  Wrel3 = Wrr23';
  Wro23=null(Wrel3); 
  
 Wteu13=null(HR13); 
 Wte13=Wteu13';
 
 Wteu23=null(HR2R13); 
 Wte23=Wteu23';
% Null at the destination Rx 
Wdess3=HR2D13*Wte23';
 Wred13=null(Wdess3'); 
 
Wdees3=HR1D23*Wte13';
 Wred23=null(Wdees3');          
%% Prpoposed Capacity -(S-R) M=2; Nr=4; Mr=6; Nd=4;       
 %% Prpoposed Scheme -(S1-R1) - Odd time slot-i 
            UU=(SNR/M *HS1R1*HS1R1'+eye(Nr));           
         CS1R1(kk)=real(log2(det(eye(Nr)+(UU))));  
%  Prpoposed Scheme of-(Ri-Di) 
          PP=SNR/Mr*S11'*S11;
        CR1D1(kk)=real(log2(det(eye(Mr)+(PP))));
        
        %% Prpoposed Scheme -(S1-R1)j
            AA=(SNR/M * Wro2*HS2R2'*HS2R2*Wro2'+Wro2*Wro2');          
         CS2R2(kk)=real(log2(det(eye(Nr)+(AA))));  
         
         %% % Prpoposed Scheme -(S1-R1) - even time slot-i 
 %% %% Prpoposed Scheme -(S1-R1) - even time slot-i 
        % %  Prpoposed Scheme of-(S-R) 
         DP=SNR/Mr* HS1R1*HS1R1';
        CS1R1EV(kk)=real(log2(det(eye(Nr)+(DP))));
%  Prpoposed Scheme of-(Ri-Di) 
        PW=SNR/Mr*Wred1'*HR1D1*Wte1'*Wte1*HR1D1'*Wred1+(Wred1'*Wred1);
        CR1D1EV(kk)=real(log2(det(eye(M)+(PW))));
                 
        %%  Prpoposed Capacity -(S-R) M1=3; Nr1=6; Mr1=9; Nd=6;
        % Prpoposed Scheme -(S1-R1) - Odd time slot-i 
%% Prpoposed Scheme -(S1-R1) - Odd time slot-i 
            UU22=(SNR/M1 *HS1R12*HS1R12'+eye(Nr1));           
         CS1R12(kk)=real(log2(det(eye(Nr1)+(UU22))));  
%  Prpoposed Scheme of-(Ri-Di) 
          PP22=SNR/Mr1*S22'*S22;
        CR1D12(kk)=real(log2(det(eye(Mr1)+(PP22))));
        
        %% Prpoposed Scheme -(S1-R1)j
            AA22=(SNR/M1 * Wro22*HS2R22'*HS2R22*Wro22'+Wro22*Wro22');          
         CS2R22(kk)=real(log2(det(eye(Nr1)+(AA22))));  
         
         %% % Prpoposed Scheme -(S1-R1) - even time slot-i 
 %% %% Prpoposed Scheme -(S1-R1) - even time slot-i 
        % %  Prpoposed Scheme of-(S-R) 
         DP22=SNR/Mr1* HS1R12*HS1R12';
        CS1R1EV2(kk)=real(log2(det(eye(Nr1)+(DP22))));
%  Prpoposed Scheme of-(Ri-Di) 
        PW22=SNR/Mr1*Wred12'*HR1D12*Wte12'*Wte12*HR1D12'*Wred12+(Wred12'*Wred12);
        CR1D1EV2(kk)=real(log2(det(eye(M1)+(PW22))));
                                      
   %% Prpoposed Capacity -(S-R)  M2=4; Nr2=8; Mr2=12; Nd=8;
   %% Prpoposed Scheme -(S1-R1) - Odd time slot-i 
           %% Prpoposed Scheme -(S1-R1) - Odd time slot-i 
            UU33=(SNR/M2 *HS1R13*HS1R13'+eye(Nr2));           
         CS1R13(kk)=real(log2(det(eye(Nr2)+(UU33))));  
%  Prpoposed Scheme of-(Ri-Di) 
          PP33=SNR/Mr2*S33'*S33;
        CR1D13(kk)=real(log2(det(eye(Mr2)+(PP33))));
        
        %% Prpoposed Scheme -(S1-R1)j
            AA33=(SNR/M2 * Wro23*HS2R23'*HS2R23*Wro23'+Wro23*Wro23');          
         CS2R23(kk)=real(log2(det(eye(Nr2)+(AA33))));  
         
         %% % Prpoposed Scheme -(S1-R1) - even time slot-i 
 %% %% Prpoposed Scheme -(S1-R1) - even time slot-i 
        % %  Prpoposed Scheme of-(S-R) 
         DP33=SNR/Mr2* HS1R13*HS1R13';
        CS1R1EV3(kk)=real(log2(det(eye(Nr2)+(DP33))));
%  Prpoposed Scheme of-(Ri-Di) 
        PW33=SNR/Mr2*Wred13'*HR1D13*Wte13'*Wte13*HR1D13'*Wred13+(Wred13'*Wred13);
        CR1D1EV3(kk)=real(log2(det(eye(M2)+(PW33))));

          end
          
     % min
       % min
        CAPS1RD1=0.5*min(CS1R1,CR1D1); %Proposed Scheme -odd time slot i
        CAPSRD12=0.5*min(CS1R1EV,CR1D1EV); %Proposed Scheme -even time slot i
        CASRDFD=(CAPS1RD1+CAPSRD12);  % M=2; Nr=6; Mr=2; Nd=4;
        
         CASRDHD= 0.5*min(CS2R2,CR1D1EV); %% HD of proposed scheme  % M=2; Nr=4; Mr=6; Nd=4;
                                           
     CAPS1RD12=0.5*min(CS1R12,CR1D12); %Proposed Scheme -odd time slot i %% M=3; Nr=6; Mr=9; Nd=6;
        CAPSRD122=0.5*min(CS1R1EV2,CR1D1EV2); %Proposed Scheme -even time slot i
        CASRDFD2=(CAPS1RD12+CAPSRD122);  % M=3; Nr=6; Mr=9; Nd=6;
        
         CASRDHD2= 0.5*min(CS2R22,CR1D1EV2); %% HD of proposed scheme    % M=2; Nr=6; Mr=2; Nd=4;
                                           
         
        CAPS1RD13=0.5*min(CS1R13,CR1D13); %Proposed Scheme -odd time slot i      %%   % M2=4; Nr2=8; Mr2=12; Nd2=8;
        CAPSRD123=0.5*min(CS1R1EV3,CR1D1EV3); %Proposed Scheme -even time slot i
        CASRDFD3=(CAPS1RD13+CAPSRD123);  % M=2; Nr=6; Mr=2; Nd=4;
        
        CASRDHD3= 0.5*min(CS2R23,CR1D1EV3); %% HD of proposed scheme  % M=2; Nr=4; Mr=6; Nd=4;
   
  
          
    % average over all channel realizations for a given SNR value
    CAPS1R1D1(SNRidx)  = mean(CASRDFD);  %% M=2; Nr=4; Mr=6; Nd=4;
    CAPSRDHDPR(SNRidx)  = mean(CASRDHD);
     
     CAPS1R1D12(SNRidx)  = mean(CASRDFD2); %% M1=3; Nr1=9; Mr1=3; Nd=6;
    CAPSRDHDPR2(SNRidx)  = mean(CASRDHD2);
     
                                        
   CAPS1R1D13(SNRidx)  = mean(CASRDFD3);  % M2=4; Nr2=12; Mr2=4; Nd=8;
   CAPSRDHDPR3(SNRidx)  = mean(CASRDHD3);
end

% plot
figure(1)

plot(SNRdBvalues, (CAPS1R1D1+CAPSRDHDPR),'k-+','linewidth',2.3)
 hold on
 plot(SNRdBvalues,(CAPS1R1D12+CAPSRDHDPR2),'k-s','linewidth',2.3)
 hold on
plot(SNRdBvalues, (CAPS1R1D13+CAPSRDHDPR3),'k:+','linewidth',2.3)

title('Sum Ergodic Capacities of HD and FD')
xlabel(' SNR(dB)')
ylabel('Sum Ergodic Capacity (bits/s/Hz)')
legend('Ms=2; Nr=4; Mr=6 ;Nd=4','Ms=3; Nr=6; Mr=9 ;Nd=6 ','Ms=4; Nr=8; Mr=12 ;Nd=8')
hold off
