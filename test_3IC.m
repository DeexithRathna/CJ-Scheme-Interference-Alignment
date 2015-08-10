clc;
close all;
clear all;

Signal_PowerdB = 0:2:20;
Signal_Powerlin = 10.^(0.1*Signal_PowerdB);
Noise_Power = 1;
Fading_Power = 1;

SNR1 = [];

% for i = 1:length(Signal_PowerdB)
    S = Signal_Powerlin(1);
%     error = 0;
%     while(error<100)
        x1 = sqrt(S) * randsrc();
        x2 = sqrt(S) * randsrc();
        x3 = sqrt(S) * randsrc();
        
        % Generation of Channel %
        H11 = randn(2,2) + randn(2,2)*1i;
        H11 = H11/norm(H11);
        
        H12 = randn(2,2) + randn(2,2)*1i;
        H12 = H12/norm(H12);
        
        H13 = randn(2,2) + randn(2,2)*1i;
        H13 = H13/norm(H13);
        
        H21 = randn(2,2) + randn(2,2)*1i;
        H21 = H21/norm(H21);
        
        H22 = randn(2,2) + randn(2,2)*1i;
        H22 = H22/norm(H11);
        
        H23 = randn(2,2) + randn(2,2)*1i;
        H23 = H23/norm(H23);
        
        H31 = randn(2,2) + randn(2,2)*1i;
        H31 = H31/norm(H31);
        
        H32 = randn(2,2) + randn(2,2)*1i;
        H32 = H32/norm(H32);
        
        H33 = randn(2,2) + randn(2,2)*1i;
        H33 = H33/norm(H33);
        
        %% Calculating Precoder Vectors
        
        [V3, D] = eig((H23^-1)*H21*(H31^-1)*H32*(H12^-1)*H13);
        V1 = (H21^-1)*H23*V3;
        V2 = (H12^-1)*H13*V3;
        
        %% Calculation of Decoder Vectors
        
        [u1 ,du1]=eig([H12*V2 H13*V3]);
        [~ ,idx]=sort(abs(diag(du1)),'ascend');
        U1=u1(:,idx(1:2));
        
        [u2 ,du2]=eig([H12*V2 H13*V3]);
        [~ ,idx]=sort(abs(diag(du2)),'ascend');
        U2=u2(:,idx(1:2));
        
        [u3 ,du3]=eig([H12*V2 H13*V3]);
        [aux ,idx]=sort(abs(diag(du3)),'ascend');
        U3=u3(:,idx(1:2));
        
        %% Received signals and Decoding
        y1 = (conj(U1))' *(H11*V1*x1 + H12*V2*x2 + H13*V3*x3);
        y2 = (conj(U2))' *(H21*V1*x1 + H22*V2*x2 + H23*V3*x3);
        y3 = (conj(U3))' *(H31*V1*x1 + H32*V2*x2 + H33*V3*x3);
        
        
        
        
%     end
% end   
