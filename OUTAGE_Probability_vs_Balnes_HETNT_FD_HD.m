clearvars;
close all;
clc;
snrdb= [-5:2:15]; %R=3; [-5:2:15]  %% with semilogy [-5:2:15] R=2
%Transmit power
% Es=2;
M=2;  % no. of antenas at source / destination 
Nr=4;  % receive antenna at relay
Mr=6;  % transmit antenna at the relay
Nd=4;  % transmit antenna at the relay
% Selected transmission rate: 2, 4 and 6 bits per transmission
% Rvalues = 4;
R=2;  %6
%bandwidth
for snri=1:length(snrdb)
    
    %no of sample for monte carlo
    count=100;   
    
    %convert snr in db to watt
    SNR=10^(snrdb(snri)/10);
    
    %initialize error counter
    er_count_raylegh=0;
    er_count_raylegh1=0;
    er_count_raylegh2=0;
    er_count_raylegh3=0;
    er_count_raylegh4=0;
    er_count_raylegh5=0;
    er_count_raylegh6=0;
    
    for i=1:count % monte carlo loop
       % generate channel realization
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
 % Null at the relay tx 
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

 %% Prpoposed Scheme -(S1-R1)i
% Prpoposed Scheme -(S1-R1) - Odd time slot-i 
            UU=(SNR/M *HS1R1*HS1R1'+eye(Nr));          
         CS1R1(i)=real(log2(det(eye(Nr)+(UU))));  
%  Prpoposed Scheme of-(Ri-Di) 
          PP=SNR/Mr*S1'*S1;
        CR1D1(i)=real(log2(det(eye(Mr)+(PP))));
                
       % Prpoposed Scheme -(S1-R1)j
            AA=(SNR/M *Wro2*HS2R2'*HS2R2*Wro2'+Wro2*Wro2');        
         CS2R2(i)=real(log2(det(eye(Nr)+(AA))));  
         
        %% %% Prpoposed Scheme -(S1-R1) - even time slot-i 
        % %  Prpoposed Scheme of-(S-R) 
         DP=SNR/Mr* HS1R1*HS1R1';
        CS1R1EV(i)=real(log2(det(eye(Nr)+(DP))));
%  Prpoposed Scheme of-(Ri-Di) 
        PW=SNR/Mr*Wred1'*HR1D1*Wte1'*Wte1*HR1D1'*Wred1+(Wred1'*Wred1);
        CR1D1EV(i)=real(log2(det(eye(M)+(PW))));
     
        %% min 
        % min
        CAPS1RD1=0.5*min(CS1R1,CR1D1); %Proposed Scheme -odd time slot i
        CAPSRD12=0.5*min(CS1R1EV,CR1D1EV); %Proposed Scheme -even time slot i
        
        CAPSRDFD=(CAPS1RD1+CAPSRD12); 
          
        %detect outage event
        if CAPSRDFD(i)<R
            er_count_raylegh=er_count_raylegh+1;
        end

 
    end

    %compute the outage probability (total error/total channel sample
%     Cnullspace(snri)=er_count_raylegh/count;
    CAPS1R1D1(snri)=er_count_raylegh/count; 
   
end

figure;

semilogy(snrdb, CAPS1R1D1,'k-+','linewidth',2.3,'MarkerSize',10)

title('Outage probability of SR-hop and RD-hops');
xlabel('SNR  (dB)');
ylabel('Outage probability')
legend('Proposed Scheme FD')