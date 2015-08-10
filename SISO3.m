%%------------------------------------------------------------------------%
%%----------#3 User SISO --- Comparision of differenct IA schemes---------%
%---------Interference Leakage Minimization Algorithm (C-J Scheme)--------%
%%------------------------------------------------------------------------%

clc;
close all;
clear all;

SNR_dB = 0:1:20;
SNR_lin = 10.^(0.1*SNR_dB);

BER1_CJ = zeros(size(SNR_dB));
BER2_CJ = BER1_CJ;
BER3_CJ = BER1_CJ;
BER_TDMA = BER1_CJ;

BER1_KT = BER1_CJ;
BER2_KT = BER1_CJ;
BER3_KT = BER1_CJ;

C1_CJ = BER1_CJ;
C2_CJ = BER2_CJ;
C3_CJ = BER3_CJ;
C1_TDMA = BER3_CJ;

N = 3;  % No of Symbol extensions = 3;
ERROR_MAX = 100;
ITER_MAX = 5000;
for i = 1:length(SNR_dB)
    error1_CJ = 0; error1_KT = 0;
    error2_CJ = 0; error2_KT = 0;
    error3_CJ = 0; error3_KT = 0;
    error1_TDMA = 0;c1_TDMA  = 0;
    cout_TDMA = 0;iter = 0;
    tic
    while( iter < ITER_MAX )
        iter = iter + 1;
        %-----------------------------------------------------------------%
        %---------------------------S-C-H-E-M-E---------------------------%
        %-----------------------------------------------------------------%
        
        %-------------Generation of channel Coefficients------------------%
        
        nT = [1,1,1];nR = [1,1,1];
        H = (1/sqrt(2))*sqrt(randn(N).^2 + randn(N).^2);
        
        %-------------No of symbol Extensions  N = 3----------------------%
        
        H11= diag(H(1,:));H12 = diag(H(2,:));H13 = diag(H(3,:));
        H = (1/sqrt(2))*sqrt(randn(N).^2 + randn(N).^2);
        H21= diag(H(1,:));H22 = diag(H(2,:));H23 = diag(H(3,:));
        H = (1/sqrt(2))*sqrt(randn(N).^2 + randn(N).^2);
        H31= diag(H(1,:));H32 = diag(H(2,:));H33 = diag(H(3,:));
        
        %------------Symbols to be transmitted----------------------------%
        
        x1 = randsrc(2,1); x2 = randsrc(); x3 = randsrc();
        X1 = sqrt(SNR_lin(i)/6)*x1;
        X2 = sqrt(SNR_lin(i)/3)*x2;
        X3 = sqrt(SNR_lin(i)/3)*x3;
        
        %------------Design of Precoding Vectors--------------------------%
        
        T = (H12)*(H21^-1)*(H23)*(H32^-1)*(H31)*(H13^-1);
        w = ones(N,1);
        A = [w,T*w];
        B = T*w;
        C = w;
        
        %------------------C-J--------S-C-H-E-M-E-------------------------%
        
        V1 = A;
        V2 = (H32^-1)*H31*C;
        V3 = (H23^-1)*H21*B;
        
        %--------------------AWGN NOISE-----------------------------------%
        
        Z1 = randn(3,1);
        Z2 = randn(3,1);
        Z3 = randn(3,1);
        
        %---------Received Signals using IA Technique---------------------%
        
        Y1 = H11*V1*X1 + H12*V2*X2 + H13*V3*X3 + Z1;
        Y2 = H21*V1*X1 + H22*V2*X2 + H23*V3*X3 + Z2;
        Y3 = H31*V1*X1 + H32*V2*X2 + H33*V3*X3 + Z3;
        
        %------------Zero Forcing Decoders Design-------------------------%
        
        U1 = [H11*V1,H12*V2,H13*V3];
        U2 = [H21*V1,H22*V2,H23*V3];
        U3 = [H31*V1,H32*V2,H33*V3];
        
        U1 = ((conj(U1))'*U1 + 0.0001*eye(size(U1'*U1)))^-1 * (conj(U1))';
        U2 = ((conj(U2))'*U2 + 0.0001*eye(size(U2'*U2)))^-1 * (conj(U2))';
        U3 = ((conj(U3))'*U3 + 0.0001*eye(size(U3'*U3)))^-1 * (conj(U3))';
        
        %--------------Symbol Recovery------------------------------------%
        
        D1 = phase(U1*Y1);
        D2 = phase(U2*Y2);
        D3 = phase(U3*Y3);
        
        D1(D1>=0) = 1;D1(D1<0) = -1;
        D2(D2>=0) = 1;D2(D2<0) = -1;
        D3(D3>=0) = 1;D3(D3<0) = -1;
        
        X1_detected = D1;
        X2_detected = D2;
        X3_detected = D3;
        
        error1_CJ = error1_CJ + sum(x1 ~= X1_detected(1:2));
        error2_CJ = error2_CJ + sum(x2 ~= X2_detected(3));
        error3_CJ = error3_CJ + sum(x3 ~= X3_detected(4));
        
        %-------------------TDMA---------SCHEME---------------------------%
        
        x1_TDMA = randsrc();
        X1_TDMA = sqrt(SNR_lin(i))*x1_TDMA;
        Y1_TDMA = H(1)*X1_TDMA + randn();
        
        %----------------Zero Forcing in TDMA-----------------------------%
        
        X1_TDMA_detected = (1/H(1))*(Y1_TDMA);
        X1_TDMA_detected(X1_TDMA_detected < 0) = -1;
        X1_TDMA_detected(X1_TDMA_detected >= 0) = 1;
        error1_TDMA = error1_TDMA + sum(X1_TDMA_detected~=x1_TDMA);
        
    end
    
    BER1_CJ(i) = error1_CJ/(2*ITER_MAX);
    BER2_CJ(i) = error2_CJ/ITER_MAX;
    BER3_CJ(i) = error3_CJ/ITER_MAX;
    BER_TDMA(i) = error1_TDMA/ITER_MAX;
    
    toc
end

%--------------------BER vs SNR dB---------------------------------%

figure;
semilogy(SNR_dB,BER1_CJ,'r--s','LineWidth',2); hold on;
semilogy(SNR_dB,BER2_CJ,'g--*','LineWidth',2); hold on;
semilogy(SNR_dB,BER3_CJ,'y--.','LineWidth',2); hold on;
semilogy(SNR_dB,BER_TDMA,'c--o','LineWidth',2);hold on;
xlabel('SNR-------> (dB)','FontSize',12,...
       'FontWeight','bold');
ylabel('BER-------->','FontSize',12,...
       'FontWeight','bold');
title('BER vs SNR 3 USER SISO (1x1)','FontSize',12,...
       'FontWeight','bold');
legend('User #1','User #2','User #3','User #1 TDMA');

DoF_CJ_ZF = (3 - BER1_CJ - BER2_CJ - BER3_CJ);
DoF_TDMA = (1-BER_TDMA);
DoF_CJ_Ideal = (3/2)* ones(size(SNR_dB));
DoF_TDMA_Ideal = ones(size(SNR_dB));

%--------------------- DoF vs SNR dB--------------------------------------%

figure;
plot(SNR_dB,DoF_CJ_ZF,'r--s','LineWidth',2); hold on;
plot(SNR_dB,DoF_TDMA,'c--o','LineWidth',2);
plot(SNR_dB,DoF_CJ_Ideal,'--p','LineWidth',2);
plot(SNR_dB,DoF_TDMA_Ideal,'--.','LineWidth',2);
xlabel('SNR------->(dB)','FontSize',12,...
       'FontWeight','bold');
ylabel('DoF------->(bits/channel Usage)','FontSize',12,...
       'FontWeight','bold');
title('DoF vs SNR 3 USER SISO (1x1)','FontSize',12,...
       'FontWeight','bold');
legend('Total (IA-CJ)','TDMA','CJ_Ideal','TDMA Ideal');

SumCapacity_CJ_ZF = DoF_CJ_ZF.*(log2(1+SNR_lin/3));
SumCapacity_TDMA = DoF_TDMA.*(log2(1+SNR_lin));
SumCapacity_TDMA_Ideal = DoF_TDMA_Ideal.*(log2(1+SNR_lin));
SumCapacity_CJ_Ideal = DoF_CJ_Ideal.*(log2(1+SNR_lin/3));

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
title('Normalized Sum Capacity vs SNR 3 USER SISO (1x1)','FontSize',12,...
       'FontWeight','bold');
legend('CJ_ZF','TDMA','CJ Ideal','TDMA Ideal');


