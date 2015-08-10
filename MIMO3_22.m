%-------------------------------------------------------------------------%
%----------------------#3 User 2x2 MIMO IA--------------------------------%
%---------Interference Leakage Minimization Algorithm (C-J Scheme)--------%
%-------------------------------------------------------------------------%

clc;
close all;
clear all;


SNR_dB = 0:1:20;
SNR_lin = 10.^(0.1*SNR_dB);

BER1_CJ = zeros(size(SNR_dB));
BER2_CJ = BER1_CJ;
BER3_CJ = BER1_CJ;
BER_TDMA = BER1_CJ;

BER1_M = BER1_CJ;
BER2_M = BER1_CJ;
BER3_M = BER1_CJ;

ERROR_MAX = 100;
ITER_MAX = 5000;

nT = [2,2,2];
nR = [2,2,2];

for i = 1:length(SNR_dB)
     error1_CJ = 0;error1_M = 0;
     error2_CJ = 0;error2_M = 0;
     error3_CJ = 0;error3_M = 0;   
     error1_TDMA = 0;
     iter = 0;
     tic
     while(iter < ITER_MAX)
         iter = iter + 1;
         
         %-----------------Channel Coefficients---------------------------%
         
         H11 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i); 
         H12 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         H13 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         
         H21 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         H22 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         H23 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         
         H31 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         H32 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         H33 = (1/sqrt(2))*(randn(2,2)  + randn(2,2)*1i);
         
         %-----------------Symbols to be transmitted----------------------%
         
         x1 = randsrc(2,1); x2 = randsrc(2,1); x3 = randsrc(2,1);
         X1 = (sqrt(SNR_lin(i)/2))*x1;
         X2 = (sqrt(SNR_lin(i)/2))*x2;
         X3 = (sqrt(SNR_lin(i)/2))*x3;
         
         
         %------------Design of Precoding Vectors-------------------------%
        
         E = (H31^-1)*H32*(H12^-1)*H13*(H23^-1)*H21;
         F = (H32^-1)*H31;
         G = (H23^-1)*H21;
         
         [V,D] = eig(E);
         [~,index] = max(max(abs(D)));
         
         V1 = V;
         V2 = F*V1; 
         V3 =  G*V1;
         
         %--------------------AWGN NOISE----------------------------------%
        
         Z1 = randn(2,1);
         Z2 = randn(2,1);
         Z3 = randn(2,1);
         
         %---------Received Signals using IA Technique--------------------%
        
         Y1 = H11*V1*X1 + H12*V2*X2 + H13*V3*X3 + Z1;
         Y2 = H21*V1*X1 + H22*V2*X2 + H23*V3*X3 + Z2;
         Y3 = H31*V1*X1 + H32*V2*X2 + H33*V3*X3 + Z3;
         
         %------------Zero Forcing Decoders Design------------------------%
        
         U1 = [H11*V1,H12*V2,H13*V3];
         U2 = [H21*V1,H22*V2,H23*V3];
         U3 = [H31*V1,H32*V2,H33*V3];
        
         U1 = ((conj(U1))'*U1 + 0.0001*eye(size(U1'*U1)))^-1 * (conj(U1))';
         U2 = ((conj(U2))'*U2 + 0.0001*eye(size(U2'*U2)))^-1 * (conj(U2))';
         U3 = ((conj(U3))'*U3 + 0.0001*eye(size(U3'*U3)))^-1 * (conj(U3))';
         
         %-------------Classical Linear MMSE Receiver---------------------%
            %----Treating Interference as coloured noise--------------%
            
            M1 = H11*V1; J1 = [H12*V2,H13*V3];
            M2 = H22*V2; J2 = [H21*V1,H23*V3];
            M3 = H33*V3; J3 = [H31*V1,H32*V2];
            
            U1_M = M1'*((1/SNR_lin(i))*eye(size(J1*J1'))+ J1*(conj(J1))'... 
                    +  M1*(conj(M1))');
            U2_M = M2'*((1/SNR_lin(i))*eye(size(J1*J1'))+ J2*(conj(J2))'...
                    +  M2*(conj(M2))');
            U3_M = M3'*((1/SNR_lin(i))*eye(size(J1*J1'))+ J3*(conj(J3))'...
                    +  M3*(conj(M3))');
            
         
         %--------------Symbol Recovery using Zero Forcing----------------%
        
         D1 = phase(U1*Y1);
         D2 = phase(U2*Y2);
         D3 = phase(U3*Y3);
        
         D1(D1>=0) = 1;D1(D1<0) = -1;
         D2(D2>=0) = 1;D2(D2<0) = -1;
         D3(D3>=0) = 1;D3(D3<0) = -1;
        
         X1_detected = D1;
         X2_detected = D2;
         X3_detected = D3;
        
         error1_CJ = error1_CJ + sum(x1 ~= X1_detected(1));
         error2_CJ = error2_CJ + sum(x2 ~= X2_detected(2));
         error3_CJ = error3_CJ + sum(x3 ~= X3_detected(3));
         
          %--------------Symbol Recovery using MMSE Receiver--------------%
        
         D1_M = phase(U1_M*Y1);
         D2_M = phase(U2_M*Y2);
         D3_M = phase(U3_M*Y3);
        
         D1_M(D1_M>=0) = 1;D1_M(D1_M<0) = -1;
         D2_M(D2_M>=0) = 1;D2_M(D2_M<0) = -1;
         D3_M(D3_M>=0) = 1;D3_M(D3_M<0) = -1;
        
         X1_Mdetected = D1_M;
         X2_Mdetected = D2_M;
         X3_Mdetected = D3_M;
        
         error1_M = error1_M + sum(x1 ~= X1_Mdetected);
         error2_M = error2_M + sum(x2 ~= X2_Mdetected);
         error3_M = error3_M + sum(x3 ~= X3_Mdetected);
         
         %-------------------TDMA---------SCHEME--------------------------%
        
         x1_TDMA = randsrc(2,1);
         X1_TDMA = sqrt(SNR_lin(i)/2)*x1_TDMA;
         Y1_TDMA = H11*X1_TDMA + randn(2,1);
         
         %----------------Zero Forcing in TDMA----------------------------%
         
         U1_TDMA = ((conj(H11))'*H11 + 0.0001*eye(size(H11'*H11)))^-1 *... 
             (conj(H11))';
         X1_TDMA_detected = (U1_TDMA)*(Y1_TDMA);
         X1_TDMA_detected(X1_TDMA_detected < 0) = -1;
         X1_TDMA_detected(X1_TDMA_detected >= 0) = 1;
         error1_TDMA = error1_TDMA + sum(X1_TDMA_detected~=x1_TDMA);
          
     end
     
     BER1_CJ(i) = error1_CJ/(2*ITER_MAX);
     BER2_CJ(i) = error2_CJ/(2*ITER_MAX);
     BER3_CJ(i) = error3_CJ/(2*ITER_MAX);
     
     BER1_M(i) = error1_M/(2*ITER_MAX);
     BER2_M(i) = error2_M/(2*ITER_MAX);
     BER3_M(i) = error3_M/(2*ITER_MAX);
     
     BER_TDMA(i) = error1_TDMA/(2*ITER_MAX);
     toc
    
end

%-------------------BER vs SNR(dB)----------------------------------------%

figure;
semilogy(SNR_dB,BER1_CJ,'r--s','LineWidth',2); hold on;
semilogy(SNR_dB,BER2_CJ,'g--*','LineWidth',2); hold on;
semilogy(SNR_dB,BER3_CJ,'y--.','LineWidth',2); hold on;
semilogy(SNR_dB,BER_TDMA,'c--o','LineWidth',2);hold on;
xlabel('SNR-------> (dB)','FontSize',12,...
       'FontWeight','bold');
ylabel('BER-------->P(error bits)','FontSize',12,...
       'FontWeight','bold');
title('BER vs SNR 3 USER MIMO (2x2)','FontSize',12,...
       'FontWeight','bold');
legend('User #1','User #2','User #3','User #1 TDMA');
 

DoF_CJ_ZF = 2*(3-BER1_CJ-BER2_CJ-BER3_CJ);
DoF_CJ_MMSE = 2*(3-BER1_M-BER2_M-BER3_M);
DoF_TDMA = 2*(1-BER_TDMA);
DoF_CJ_Ideal = 3*ones(size(SNR_dB));
DoF_TDMA_Ideal = 2*ones(size(SNR_dB));

%--------------Degrees of Freedom vs SNR dB-------------------------------%

figure;
plot(SNR_dB,DoF_CJ_ZF,'r--s','LineWidth',2); hold on;
plot(SNR_dB,DoF_CJ_MMSE,'g--.','LineWidth',2); hold on;
plot(SNR_dB,DoF_TDMA,'c--o','LineWidth',2);
plot(SNR_dB,DoF_CJ_Ideal,'--p','LineWidth',2);
plot(SNR_dB,DoF_TDMA_Ideal);

xlabel('SNR ------- >(dB)','FontSize',12,...
       'FontWeight','bold');
ylabel('DoF ------- >(bits/channel Usage)','FontSize',12,...
       'FontWeight','bold');
title('DoF vs SNR 3 USER MIMO (2x2)','FontSize',12,...
       'FontWeight','bold');
legend('(IA-CJ-ZF)','(CJ-MMSE)','TDMA','IA Ideal','TDMA Ideal');


SumCapacity_CJ_ZF = DoF_CJ_ZF.*(log2(1+SNR_lin));
SumCapacity_CJ_MMSE = DoF_CJ_MMSE.*(log2(1+SNR_lin));
SumCapacity_TDMA = DoF_TDMA.*(log2(1+SNR_lin));
SumCapacity_TDMA_Ideal = DoF_TDMA_Ideal.*(log2(1+SNR_lin));
SumCapacity_CJ_Ideal = DoF_CJ_Ideal.*(log2(1+SNR_lin));

%-----------------------Sum Capacity vs SNR dB----------------------------%

figure;
plot(SNR_dB,SumCapacity_CJ_ZF,'r--s','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_CJ_MMSE,'--.','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_TDMA,'c--o','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_CJ_Ideal,'--p','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_TDMA_Ideal,'--.','LineWidth',2);hold on;
xlabel('SNR dB  ------- >','FontSize',12,...
       'FontWeight','bold');
ylabel('Normalized Capacity (bits/Hz) ---------- >','FontSize',12,...
       'FontWeight','bold');
title('Normalized Sum Capacity vs SNR  3 USER MIMO (2x2)','FontSize',12,...
       'FontWeight','bold');
legend('CJ ZF','CJ MMSE','TDMA','CJ Ideal','TDMA Ideal');
   
 

