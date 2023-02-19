%%% Please cite this article "Half-Duplex and Full-Duplex Interference
%%% Mitigation in Relays Assisted Heterogeneous Network" PlosOne ID: PONE-D-22-34295
%by Moubachir Madani Fadoul and Chee-Onn Chow"
%%% Figure number 7.


clearvars;
close all;
clc;
% number of channel realization
It = 10;
snrdb=2; %2;

% Number of transmit and receive antennas (2-4-6-6)or (2-6-2-4) or
% (2-2-4-6)- %%% HERE antt(2-4-6-4)
M=2;  % no. of antenas at source / destination 
Nr=4;  % receive antenna at relay
Mr=6;  % transmit antenna at the relay
Nd=4;  % transmit antenna at the relay

% Selected transmission rate: 2, 4 and 6 bits per transmission
R=[5:1:8];    %[2:2:8]   R=[0:2:10]; [6:1:11]
RXX=[5:1:8];   %[2:2:8]  RXX=[0:1:5];[5:1:10];  
% RYY=[0:2:10];  
% RZZ=[0:2:10]; 
RAA=[5:1:8];  %RAA=[0:2:10]; [6:1:11];
RBB=[5:1:8];    %RBB=[0:1:5]; [5:1:10];
RCC=[5:1:8];

R0=1;  R1=2;  R2=3;  R3=4;   %R0=0; R1=2; R2=4; R3=6; R4=8; R5=10;

 %convert snr in db to watt
    snr=10^(snrdb/10);

%  for snri=1:length(snr)
    
    %no of sample for monte carlo
    count=100;   

    %initialize error counter
    er_count_raylegh0=0; er_count_raylegh1=0; er_count_raylegh2=0;  er_count_raylegh3=0; 
    
    SNRidx=0;
     for i=1:count % monte carlo loop
        SNRidx=SNRidx+1;
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
  Wtu1=null(HR1); 
  Wto1=Wtu1';
  
  Wrr2=HR1R2*Wto1';
  Wrel = Wrr2';
  Wro2=null(Wrel); 
  
%  PSHRR2 =pinv(HR1);  
% Wtu1=(PSHRR2*HR1R2*Wt2');  
% Wt1=Wtu1';

 Wteu1=null(HR1); 
 Wte1=Wteu1';
 
 Wteu2=null(HR2R1); 
 Wte2=Wteu2';
% Null at the destination Rx 
Wdess=HR2D1*Wte2';
 Wred1=null(Wdess'); 
 
Wdees=HR1D2*Wte1';
Wred2=null(Wdees'); 

%% Prpoposed Scheme -(S1-R1) - Odd time slot-i 
            UU=(snr/M *HS1R1*HS1R1'+eye(Nr));
%               RR=SNR/Mr *(Wr1'*HR1'*HR1*Wr1);  %RSI ->becomes zero
%               YY=SNR/Mr *(Wr1'*HR2D1'*Wt2'*Wt2*HR2D1*Wr1');  %RDI  ->zeros            
         CS1R1(i)=real(log2(det(eye(Nr)+(UU))));  
%  Prpoposed Scheme of-(Ri-Di) 
%          PP=SNR/Mr*HR1D1'*HR1D1;
          PP=snr/Mr*S1'*S1;
        CR1D1(i)=real(log2(det(eye(Mr)+(PP))));
                
       % Prpoposed Scheme -(S1-R1)j
            AA=(snr/M *Wro2*HS2R2'*HS2R2*Wro2'+Wro2*Wro2');
% %             DR=SNR/Mr *(HR1*Wt2'*Wt2*HR1');  %RSI ->becomes zero
% %               YY=SNR/Mr *(Wr1'*HR2D1'*Wt2'*Wt2*HR2D1*Wr1');  %RDI  ->zeros          
         CS2R2(i)=real(log2(det(eye(Nr)+(AA))));  
         
        %% %% Prpoposed Scheme -(S1-R1) - even time slot-i 
        % %  Prpoposed Scheme of-(S-R) 
         DP=snr/Mr* HS1R1*HS1R1';
        CS1R1EV(i)=real(log2(det(eye(Nr)+(DP))));
%  Prpoposed Scheme of-(Ri-Di) 
%          PP=SNR/Mr*HR1D1'*HR1D1;
%           PW=SNR/Mr*Wred1'*HR1D1'*HR1D1*Wred1;
        PW=snr/Mr*Wred1'*HR1D1*Wte1'*Wte1*HR1D1'*Wred1+(Wred1'*Wred1);
        CR1D1EV(i)=real(log2(det(eye(M)+(PW))));
      
%               PWxx=SNR/Mr*Wred2'*HR2D2*Wte2'*Wte2*HR2D2'*Wred2+(Wred2'*Wred2);
%         CR2D2EV(kk)=real(log2(det(eye(M)+(PWxx))));
  
                
        %% 
     
   % min
        CAPS1RD1=0.5*min(CS1R1,CR1D1); %Proposed Scheme -odd time slot i
        CAPSRD12=0.5*min(CS1R1EV,CR1D1EV); %Proposed Scheme -even time slot i
        CAPSRDFD=(CAPS1RD1+CAPSRD12); 
          
       CASRDHD= 0.5*(CS2R2+CR1D1EV); %% HD of proposed scheme                           

% % % % %               %detect outage event
                  if CAPSRDFD(i)<R0
             er_count_raylegh0=er_count_raylegh0+1;
                 end
                  if CAPSRDFD(i)<R1
             er_count_raylegh1=er_count_raylegh1+1;
                  end
                  if CAPSRDFD(i)<R2
             er_count_raylegh2=er_count_raylegh2+1;
                  end
                  if CAPSRDFD(i)<R3
             er_count_raylegh3=er_count_raylegh3+1;
                  end
               
                            
end                       

% %compute the outage probability (total error/total channel sample
    CSVDNU0=er_count_raylegh0/count;
    CSVDNU1=er_count_raylegh1/count;
    CSVDNU2=er_count_raylegh2/count;
    CSVDNU3=er_count_raylegh3/count;
   
 CSVDNNALL = [CSVDNU0,CSVDNU1,CSVDNU2,CSVDNU3];

% plot figure;
semilogy(R,CSVDNNALL,'k:s','linewidth',2.3, 'MarkerSize',10);  %%k-s
 
title('Outage probability of SR-hop and RD-hop');
xlabel('Rate (Bits/s/Hz)');
ylabel('Outage Probability')
legend('Proposed Scheme HD')
hold off