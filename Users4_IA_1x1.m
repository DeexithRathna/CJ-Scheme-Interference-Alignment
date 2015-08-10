%-------------------------------------------------------------------------%
%----------------------#4 User 1x1 SISO IA--------------------------------%
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
BER4_CJ = BER1_CJ;
BER_TDMA = BER1_CJ;

ERROR_MAX = 100;
ITER_MAX = 5000;

nT = [1,1,1];          % ----- # Tx Antenna at each user ---------------- %
nR = [1,1,1];          % ----- # Rx Antenna at each user ---------------- % 
n = 1;
K = 4;                 % ------ No of users------------------------------ %
N = (K-1)*(K-2) - 1;   
Mn = (n+1)^N + n^N;    % ---------No of Signal Extensions --------------- %
w = ones(33,1);

for i = 1:length(SNR_dB)
     error1_CJ = 0;
     error2_CJ = 0;
     error3_CJ = 0;
     error4_CJ = 0;
     error1_TDMA = 0;
     iter = 0;
     tic
     while(iter < ITER_MAX)
         iter = iter + 1;
         
         %-----------------Channel Coefficients---------------------------%
         
         H = (1/sqrt(2))*((randn(33,K^2) + randn(33,K^2)*1i));
         H11 = diag(H(:,1)); H12 = diag(H(:,2)); 
         H13 = diag(H(:,3)); H14 = diag(H(:,4));
         H21 = diag(H(:,5)); H22 = diag(H(:,6));
         H23 = diag(H(:,7)); H24 = diag(H(:,8));
         H31 = diag(H(:,9)); H32 = diag(H(:,2));
         H33 = diag(H(:,11));H34 = diag(H(:,12));
         H41 = diag(H(:,13));H42 = diag(H(:,14));
         H43 = diag(H(:,15));H44 = diag(H(:,16));
         
         
         %-----------------Symbols to be transmitted----------------------%
         
         x1 = randsrc(32,1); x2 = randsrc(); x3 = randsrc(); x4 = randsrc();
         X1 = sqrt(SNR_lin(i)/64)*x1;
         X2 = sqrt(SNR_lin(i)/32)*x2;
         X3 = sqrt(SNR_lin(i)/32)*x3;
         X4 = sqrt(SNR_lin(i)/32)*x4;
         
         %------------Design of Precoding Vectors-------------------------%
         
         S2 = (H12^-1)*H13*(H23^-1)*H21;
         S3 = (H13^-1)*H13*(H23^-1)*H21;
         S4 = (H14^-1)*H13*(H23^-1)*H21;
         T23 = (H21^-1)*H23*S3;
         T24 = (H21^-1)*H24*S4;
         T32 = (H31^-1)*H32*S2;
         T34 = (H31^-1)*H34*S4;
         T42 = (H41^-1)*H42*S2;
         T43 = (H41^-1)*H43*S3;
         B =  w;
         V1 = [w,T24*w,T32*w,T34*w,T42*w,T43*w,...
                 T24*T32*w,T24*T34*w,T24*T42*w,T24*T43*w,...
                 T32*T34*w,T32*T42*w,T32*T43*w,...
                 T34*T42*w,T34*T43*w, T42*T43*w...
                 T24*T32*T34*w,T24*T32*T42*w,T24*T32*T43*w,...
                 T24*T34*T42*w, T24*T34*T43*w, T24*T42*T43*w,...
                 T32*T34*T42*w,T32*T34*T43*w,T32*T42*T43*w,T34*T42*T43*w,...
                 T24*T32*T34*T42*w,T24*T32*T34*T43*w,...
                 T24*T32*T42*T43*w,T24*T34*T42*T43*w,...
                 T32*T34*T42*T43*w,T24*T32*T34*T42*T43*w];
         V2 = S2*B;
         V3 = S3*B;
         V4 = S4*B;
         
         %--------------------AWGN NOISE----------------------------------%
         
         Z1 = (1/sqrt(2))*(randn(33,1) + randn(33,1)*1i);
         Z2 = (1/sqrt(2))*(randn(33,1) + randn(33,1)*1i);
         Z3 = (1/sqrt(2))*(randn(33,1) + randn(33,1)*1i);
         Z4 = (1/sqrt(2))*(randn(33,1) + randn(33,1)*1i);
        
         %---------Received Signals using IA Technique--------------------%
        
         Y1 = H11*V1*X1 + H12*V2*X2 + H13*V3*X3 + H14*V4*X4 + Z1;
         Y2 = H21*V1*X1 + H22*V2*X2 + H23*V3*X3 + H24*V4*X4 + Z2;
         Y3 = H31*V1*X1 + H32*V2*X2 + H33*V3*X3 + H34*V4*X4 + Z3;
         Y4 = H41*V1*X1 + H42*V2*X2 + H43*V3*X3 + H44*V4*X4 + Z4;
         
         %------------Zero Forcing Decoders Design------------------------%
        
         U1 = [H11*V1,H12*V2,H13*V3, H14*V4];
         U2 = [H21*V1,H22*V2,H23*V3, H24*V4];
         U3 = [H31*V1,H32*V2,H33*V3, H34*V4];
         U4 = [H41*V1,H42*V2,H43*V3, H44*V4];
        
         U1 = ((conj(U1))'*U1 + 0.01*eye(size(U1'*U1)))^-1 * (conj(U1))';
         U2 = ((conj(U2))'*U2 + 0.01*eye(size(U2'*U2)))^-1 * (conj(U2))';
         U3 = ((conj(U3))'*U3 + 0.01*eye(size(U3'*U3)))^-1 * (conj(U3))';
         U4 = ((conj(U4))'*U4 + 0.01*eye(size(U4'*U4)))^-1 * (conj(U4))';
         
         %--------------Symbol Recovery using Zero Forcing----------------%
        
         D1 = phase(U1*Y1);
         D2 = phase(U2*Y2);
         D3 = phase(U3*Y3);
         D4 = phase(U4*Y4);
        
         D1(D1>=0) = 1;D1(D1<0) = -1;
         D2(D2>=0) = 1;D2(D2<0) = -1;
         D3(D3>=0) = 1;D3(D3<0) = -1;
         D4(D4>=0) = 1;D4(D4<0) = -1;
        
         X1_detected = D1;
         X2_detected = D2;
         X3_detected = D3;
         X4_detected = D4;
        
         error1_CJ = error1_CJ + sum(x1 ~= X1_detected(1:32));
         error2_CJ = error2_CJ + sum(x2 ~= X2_detected(33));
         error3_CJ = error3_CJ + sum(x3 ~= X3_detected(34));
         error4_CJ = error4_CJ + sum(x4 ~= X4_detected(35));
         
         %-------------------TDMA---------SCHEME--------------------------%
        
         x1_TDMA = randsrc();
         X1_TDMA = sqrt(SNR_lin(i))*x1_TDMA;
         h = (1/sqrt(2))*(randn() + randn()*1i);
         Y1_TDMA = h*X1_TDMA + randn();
         
         %----------------Zero Forcing in TDMA----------------------------%
         
        X1_TDMA_detected = (1/h)*(Y1_TDMA);
        X1_TDMA_detected(X1_TDMA_detected < 0) = -1;
        X1_TDMA_detected(X1_TDMA_detected >= 0) = 1;
        error1_TDMA = error1_TDMA + sum(X1_TDMA_detected~=x1_TDMA);
          
     end
     
     BER1_CJ(i) = error1_CJ/(32*ITER_MAX);
     BER2_CJ(i) = error2_CJ/(ITER_MAX);
     BER3_CJ(i) = error3_CJ/(ITER_MAX);
     BER4_CJ(i) = error4_CJ/ITER_MAX;

     BER_TDMA(i) = error1_TDMA/(2*ITER_MAX);
     toc
    
end

%-------------------BER vs SNR(dB)----------------------------------------%

figure;
semilogy(SNR_dB,BER1_CJ,'r--s','LineWidth',2); hold on;
semilogy(SNR_dB,BER2_CJ,'g--*','LineWidth',2); hold on;
semilogy(SNR_dB,BER3_CJ,'y--.','LineWidth',2); hold on;
semilogy(SNR_dB,BER4_CJ); hold on;
semilogy(SNR_dB,BER_TDMA,'c--o');hold on;
xlabel('SNR-------> (dB)','FontSize',12,...
       'FontWeight','bold');
ylabel('BER-------->P(error bits)','FontSize',12,...
       'FontWeight','bold');
title('BER vs SNR 3 USER MIMO (2x2)','FontSize',12,...
       'FontWeight','bold');
legend('User #1','User #2','User #3','User #4','User #1 TDMA');

DoF_CJ_ZF = (4-BER1_CJ-BER2_CJ-BER3_CJ-BER4_CJ);
DoF_TDMA = (1-BER_TDMA);
DoF_CJ_Ideal = 2*ones(size(SNR_dB));
DoF_TDMA_Ideal = 1*ones(size(SNR_dB));

%--------------Degrees of Freedom vs SNR dB-------------------------------%

figure;
plot(SNR_dB,DoF_CJ_ZF,'r--s','LineWidth',2); hold on;
plot(SNR_dB,DoF_TDMA,'c--o','LineWidth',2);
plot(SNR_dB,DoF_CJ_Ideal,'--p','LineWidth',2);
plot(SNR_dB,DoF_TDMA_Ideal);

xlabel('SNR ------- >(dB)','FontSize',12,...
       'FontWeight','bold');
ylabel('DoF ------- >(bits/channel Usage)','FontSize',12,...
       'FontWeight','bold');
title('DoF vs SNR 4 USER SISO (1x1)','FontSize',12,...
       'FontWeight','bold');
legend('(IA-CJ-ZF)','TDMA','IA Ideal','TDMA Ideal');


SumCapacity_CJ_ZF = DoF_CJ_ZF.*(log2(1+SNR_lin));
SumCapacity_TDMA = DoF_TDMA.*(log2(1+SNR_lin));
SumCapacity_TDMA_Ideal = DoF_TDMA_Ideal.*(log2(1+SNR_lin));
SumCapacity_CJ_Ideal = DoF_CJ_Ideal.*(log2(1+SNR_lin));

%-----------------------Sum Capacity vs SNR dB----------------------------%

figure;
plot(SNR_dB,SumCapacity_CJ_ZF,'r--s','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_TDMA,'c--o','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_CJ_Ideal,'--p','LineWidth',2);hold on;
plot(SNR_dB,SumCapacity_TDMA_Ideal,'--.','LineWidth',2);hold on;
xlabel('SNR dB  ------- >','FontSize',12,...
       'FontWeight','bold');
ylabel('Normalized Capacity (bits/Hz) ---------- >','FontSize',12,...
       'FontWeight','bold');
title('Normalized Sum Capacity vs SNR  3 USER MIMO (2x2)','FontSize',12,...
       'FontWeight','bold');
legend('CJ ZF','TDMA','CJ Ideal','TDMA Ideal');

