function [DS, BU, PC, MU, GN, QN, A_finalre] = functionComputeSE_AP_uplink_anal_ADC_FZF(Psi,R,K,L,D,p_uplink,alpha,pilotIndex, p_pilot,tau_p,N,c_matrix,gamma_matrix, beta_matrix, theta_matrix)
%% Prepare to save analytical results
DS = zeros(K,2,2)  ;
BU = zeros(K,2,2)  ;
PC = zeros(K,2,2)  ;
MU = zeros(K,2,2)  ;
GN = zeros(K,2,2)  ;
QN = zeros(K,2,2)  ;
A  = zeros(L,K,2);
A_FZF = zeros(L,K,2);
A_finalre = zeros(L,K,2,2);
A_part1 = zeros(L,L,K);
A_part1_FZF = zeros(L,L,K);
A_part2 = zeros(L,L,K,K);
A_part2_FZF = zeros(L,L,K);
A_part3 = zeros(L,K);
A_part3_FZF = zeros(L,K);

for k = 1:K
    %% MR
    % Level 2
    A(:,k,1) = (1/L)*ones(L,1);
    % Level 3
    A_A = ones(L,K,2);
    for l = 1:L
        A_BU = computeBU(A_A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, pilotIndex);
        [~, A_PC] = computePC(A_A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, pilotIndex, L, K);
        A_MU = computeMU(A_A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, pilotIndex, L, K);
        A_GN = computeGN(A_A, D, R, Psi, alpha, p_pilot, k, tau_p, L, K, pilotIndex);
        A_QN = computeQN(A_A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, N, pilotIndex);
        A_part1(l,l,k) = A_BU(1,l) + A_PC(1,l) + A_MU(1,l) + A_QN(1,l) + A_GN(1,l);
    end
    
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                b_k = computeA_part2(D, R, Psi, alpha, p_pilot, k, k1, tau_p, L, pilotIndex);
                A_part2(:,:,k,k1) = alpha^2*p_uplink(k1)*b_k*b_k';
            end
        end
    end
    A_part1(:,:,k) = A_part1(:,:,k) + squeeze(sum(A_part2(:,:,k,:),4));
    A_part3(:,k) = computeA_part2(D, R, Psi, alpha, p_pilot, k, k, tau_p, L, pilotIndex);
    temp_1 = A_part1(:,:,k);
    temp_1(all(temp_1==0,2),:)=[];
    % %     temp_1(:,all(temp_1==0,1))=[];
    temp_2 = A_part3(:,k);
    temp_2(find(temp_2==0))=[];
    A(:,k,2) = sqrt(alpha^2*p_uplink(k))*(temp_1\ temp_2);
    
    %% FZF
    % Level 2
    A_FZF(:,k,1) = (1/L)*ones(L,1);
    % Level 3
    for l = 1: L
        A_BU_FZF = computeBU_FZF(alpha, p_uplink, A_A, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p, L);
        [~, A_PC_FZF] = computePC_FZF(alpha, c_matrix, pilotIndex, p_uplink, A_A, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p,K,L);
        A_MU_FZF = computeMU_FZF(alpha, pilotIndex, p_uplink, A_A, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p,K,L);
        A_GN_FZF = computeGN_FZF(alpha, A_A, N, tau_p, theta_matrix, k,K,L);
        A_QN_FZF = computeQN_FZF(alpha, A_A, N, tau_p, theta_matrix, k, beta_matrix,K,L,p_uplink);
        A_part1_FZF(l,l,k) = A_BU_FZF(1,l) + A_PC_FZF(1,l) + A_MU_FZF(1,l) + A_QN_FZF(1,l) + A_GN_FZF(1,l);
    end
    
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                b_k_vec = c_matrix(:,k1);
                A_part2_FZF(:,:,k) = A_part2_FZF(:,:,k) + alpha^2*p_uplink(k1)*b_k_vec*b_k_vec';
            end
        end
    end
    A_part1_FZF(:,:,k) = A_part1_FZF(:,:,k) + A_part2_FZF(:,:,k);
    A_part3_FZF(:,k) = c_matrix(:,k);
    temp_1 = A_part1_FZF(:,:,k);
    temp_1(all(temp_1==0,2),:)=[];
    %     temp_1(:,all(temp_1==0,1))=[];
    temp_2 = A_part3_FZF(:,k);
    temp_2(find(temp_2==0))=[];
    A_FZF(:,k,2) = sqrt(alpha^2*p_uplink(k))*(temp_1\ temp_2);
    
    %%
    DS(k,:,1) = computeDS(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, pilotIndex)        ;
    BU(k,:,1) = sum(computeBU(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, pilotIndex),2)   ;
    [PC_all, ~] = computePC(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, pilotIndex, L, K)    ;
    PC(k,:,1) = PC_all;
    MU(k,:,1) = sum(computeMU(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, pilotIndex, L, K),2)   ;
    GN(k,:,1) = sum(computeGN(A, D, R, Psi, alpha, p_pilot, k, tau_p, L, K, pilotIndex),2)             ;
    QN(k,:,1) = sum(computeQN(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, N, pilotIndex),2);
    
    DS(k,:,2)  = computeDS_FZF(alpha, p_uplink, A_FZF, c_matrix,k, L);
    
    temp_22 = computeBU_FZF(alpha, p_uplink, A_FZF, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p, L);
    BU(k,:,2) = sum(temp_22,2);
    
    [temp_23,~] = computePC_FZF(alpha, c_matrix, pilotIndex, p_uplink, A_FZF, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p,K,L);
    PC(k,:,2) = sum(temp_23,2);
    
    temp_24 = computeMU_FZF(alpha, pilotIndex, p_uplink, A_FZF, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p,K,L);
    MU(k,:,2) = sum(temp_24,2);
    
    temp_25 = computeGN_FZF(alpha, A_FZF, N, tau_p, theta_matrix, k,K,L);
    GN(k,:,2) = sum(temp_25,2);
    temp_26 = computeQN_FZF(alpha, A_FZF, N, tau_p, theta_matrix, k, beta_matrix,K,L,p_uplink);
    QN(k,:,2) = sum(temp_26,2);
    
    A_finalre(:,k,:,1) = A(:,k,:);
    A_finalre(:,k,:,2) = A_FZF(:,k,:);
    
end
end

%% sub-function
function DS = computeDS(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, pilotIndex)

DS = zeros(2,1);

for l = 1:L
    DS(1,1) = DS(1,1)+(conj(A(l,k,1)))*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
    DS(2,1) = DS(2,1)+(conj(A(l,k,2)))*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
end

DS(1,1) = alpha^6*(p_pilot*tau_p)^2*p_uplink(k)*abs(DS(1,1))^2;
DS(2,1) = alpha^6*(p_pilot*tau_p)^2*p_uplink(k)*abs(DS(2,1))^2;
end

function BU = computeBU(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, pilotIndex)

BU = zeros(2,L);

for l = 1:L
    BU(1,l) = alpha^4*p_pilot*tau_p*p_uplink(k)*abs(A(l,k,1))^2*trace(R(:,:,l,k)*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
    BU(2,l) = alpha^4*p_pilot*tau_p*p_uplink(k)*abs(A(l,k,2))^2*trace(R(:,:,l,k)*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
end

end

function [PC_all, PC_A ] = computePC(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, pilotIndex, L, K)

PC_all = zeros(2,1);
PC_part = zeros(2,K,L);
PC1= zeros(2,K,L);
PC2= zeros(2,K,L);
PC3= zeros(2,K);
PC4= zeros(2,K,L);


for i = 1:K
    if pilotIndex(i)==pilotIndex(k)
        if i~=k
            for l = 1:L
                
                PC1(1,i,l) = alpha^2*p_uplink(i)*abs(A(l,k,1))^2*(alpha^2*p_pilot*tau_p)^2*abs(trace((R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k))^0.5*D(:,:,l,k)*(R(:,:,l,i)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,i))^0.5))^2;
                PC2(1,i,l) = alpha^2*p_uplink(i)*abs(A(l,k,1))^2*(alpha^2*p_pilot*tau_p)*trace(R(:,:,l,i)*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
                PC4(1,i,l) = alpha^2*p_uplink(i)*(alpha^2*p_pilot*tau_p)^2*abs(A(l,k,1))^2*abs(trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,i)))^2;
                
                PC1(2,i,l) = alpha^2*p_uplink(i)*abs(A(l,k,2))^2*(alpha^2*p_pilot*tau_p)^2*abs(trace((R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k))^0.5*D(:,:,l,k)*(R(:,:,l,i)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,i))^0.5))^2;
                PC2(2,i,l) = alpha^2*p_uplink(i)*abs(A(l,k,2))^2*(alpha^2*p_pilot*tau_p)*trace(R(:,:,l,i)*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
                PC4(2,i,l) = alpha^2*p_uplink(i)*(alpha^2*p_pilot*tau_p)^2*abs(A(l,k,2))^2*abs(trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,i)))^2;
                
                PC3(1,i) = PC3(1,i) + (alpha^2*p_pilot*tau_p)*(conj(A(l,k,1)))*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,i));
                PC3(2,i) = PC3(2,i) + (alpha^2*p_pilot*tau_p)*(conj(A(l,k,2)))*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,i));
                
                PC_part(1,i,l) = (PC2(1,i,l));
                PC_part(2,i,l) = (PC2(2,i,l));
            end
            PC3(1,i) =  alpha^2*p_uplink(i)*abs(PC3(1,i))^2;
            PC3(2,i) =  alpha^2*p_uplink(i)*abs(PC3(1,i))^2;
        end
    end
end
PC_A = squeeze(sum(PC_part,2));

PC_all(1,1) = (squeeze(sum(sum(PC1(1,:,:))))+squeeze(sum(sum(PC2(1,:,:))))+sum(PC3(1,:))-squeeze(sum(sum(PC4(1,:,:)))));
PC_all(2,1) = (squeeze(sum(sum(PC1(2,:,:))))+squeeze(sum(sum(PC2(2,:,:))))+sum(PC3(1,:))-squeeze(sum(sum(PC4(2,:,:)))));
end

function MU = computeMU(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, pilotIndex, L, K)

MU = zeros(2,L);

for l = 1:L
    for i = 1:K
        if pilotIndex(i)~=pilotIndex(k)
            
            MU(1,l) = MU(1,l) + alpha^4*p_pilot*tau_p*p_uplink(i)*abs(A(l,k,1))^2*trace(R(:,:,l,i)*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
            MU(2,l) = MU(2,l) + alpha^4*p_pilot*tau_p*p_uplink(i)*abs(A(l,k,2))^2*trace(R(:,:,l,i)*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
        end
    end
end
end

function GN = computeGN(A, D, R, Psi, alpha, p_pilot, k, tau_p, L, K, pilotIndex)

GN = zeros(2,L);

for l = 1:L
    GN(1,l) = alpha^4*p_pilot*tau_p*abs(A(l,k,1))^2*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
    GN(2,l) = alpha^4*p_pilot*tau_p*abs(A(l,k,2))^2*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
end
end

function QN = computeQN(A, D, R, Psi, alpha, p_pilot, p_uplink, k, tau_p, L, K, N, pilotIndex)

QN = zeros(2,L);
temp = zeros(N,N,K,K);

for l = 1:L
    for k1 = 1:K
        temp(:,:,k,k1) = p_uplink(k1)*R(:,:,l,k1);
    end
    QN(1,l) = alpha^3*(1-alpha)*p_pilot*tau_p*abs(A(l,k,1))^2*trace(diag(diag(squeeze(sum(temp(:,:,k,:),4))+eye(N)))*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
    QN(2,l) = alpha^3*(1-alpha)*p_pilot*tau_p*abs(A(l,k,2))^2*trace(diag(diag(squeeze(sum(temp(:,:,k,:),4))+eye(N)))*D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k));
end
end

function A_part2 = computeA_part2(D, R, Psi, alpha, p_pilot, k, k1, tau_p, L, pilotIndex)

A_part2 = zeros(L,1);

for l = 1:L
    A_part2(l,1) = alpha^2*p_pilot*tau_p*trace(D(:,:,l,k)*R(:,:,l,k)/squeeze(Psi(:,:,l,pilotIndex(k)))*R(:,:,l,k1));
end
end

function DS_FZF = computeDS_FZF(alpha, p_uplink, A_FZF, c_matrix,k, L)

DS_FZF = zeros(2,1);
DS_FZF_temp_1 = 0;
DS_FZF_temp_2 = 0;
for l = 1:L
    DS_FZF_temp_1 = DS_FZF_temp_1 + (conj(A_FZF(l,k,1))*c_matrix(l,k));
    DS_FZF_temp_2 = DS_FZF_temp_2 + (conj(A_FZF(l,k,2))*c_matrix(l,k));
end
DS_FZF(1,1) = alpha^2*p_uplink(k)*abs(DS_FZF_temp_1)^2;
DS_FZF(2,1) = alpha^2*p_uplink(k)*abs(DS_FZF_temp_2)^2;
end

function BU_FZF = computeBU_FZF(alpha, p_uplink, A_FZF, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p, L)

BU_FZF = zeros(2,L);
for l = 1:L
        BU_FZF(1,l) = BU_FZF(1,l)+alpha^2*p_uplink(k)*conj(A_FZF(l,k,1))*conj(A_FZF(l,k,1))*(beta_matrix(l,k)-gamma_matrix(l,k))/((N-tau_p)*theta_matrix(l,k));
        BU_FZF(2,l) = BU_FZF(2,l)+alpha^2*p_uplink(k)*conj(A_FZF(l,k,2))*conj(A_FZF(l,k,2))*(beta_matrix(l,k)-gamma_matrix(l,k))/((N-tau_p)*theta_matrix(l,k));
end
end

function [PC_FZF_all,PC_FZF_part]  = computePC_FZF(alpha, c_matrix, pilotIndex, p_uplink, A_FZF, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p, K,L)

PC_FZF_part = zeros(2,L);
PC_FZF_all = zeros(2,1);
PC_FZF_temp_1_2 = 0;
PC_FZF_temp1 = 0;
PC_FZF_temp_2_2 = 0;
PC_FZF_temp2 = 0;

for k1 = 1:K
    if pilotIndex(k1)==pilotIndex(k)
        if k1~=k
            PC_FZF_temp_1_1 = 0;
            PC_FZF_temp_2_1 = 0;
            for l = 1:L
                PC_FZF_temp_1_1 = PC_FZF_temp_1_1 + (conj(A_FZF(l,k,1))*c_matrix(l,k1));
                PC_FZF_temp_2_1 = PC_FZF_temp_2_1 + (conj(A_FZF(l,k,2))*c_matrix(l,k1));
                PC_FZF_temp_1_2 = PC_FZF_temp_1_2 + alpha^2*p_uplink(k1)*conj(A_FZF(l,k,1))*conj(A_FZF(l,k,1))*(beta_matrix(l,k1)-gamma_matrix(l,k1))/((N-tau_p)*theta_matrix(l,k));
                PC_FZF_temp_2_2 = PC_FZF_temp_2_2 + alpha^2*p_uplink(k1)*conj(A_FZF(l,k,2))*conj(A_FZF(l,k,2))*(beta_matrix(l,k1)-gamma_matrix(l,k1))/((N-tau_p)*theta_matrix(l,k));
            end
            PC_FZF_temp1 = PC_FZF_temp1 + alpha^2*p_uplink(k1)*(abs(PC_FZF_temp_1_1)^2);
            PC_FZF_temp2 = PC_FZF_temp2 + alpha^2*p_uplink(k1)*(abs(PC_FZF_temp_2_1)^2);
        end
    end
end


for k1 = 1:K
    if pilotIndex(k1)==pilotIndex(k)
        if k1~=k
            for l = 1:L
                PC_FZF_part(1,l) = PC_FZF_part(1,l)+alpha^2*p_uplink(k1)*conj(A_FZF(l,k,1))*conj(A_FZF(l,k,1))*(beta_matrix(l,k1)-gamma_matrix(l,k1))/((N-tau_p)*theta_matrix(l,k));
                PC_FZF_part(2,l) = PC_FZF_part(2,l)+alpha^2*p_uplink(k1)*conj(A_FZF(l,k,2))*conj(A_FZF(l,k,2))*(beta_matrix(l,k1)-gamma_matrix(l,k1))/((N-tau_p)*theta_matrix(l,k));
            end
        end
    end
end
PC_FZF_all(1,1) = (PC_FZF_temp1 + PC_FZF_temp_1_2);
PC_FZF_all(2,1) = (PC_FZF_temp2 + PC_FZF_temp_2_2);
end


function MU_FZF = computeMU_FZF(alpha, pilotIndex, p_uplink, A_FZF, beta_matrix, gamma_matrix, theta_matrix, k, N, tau_p,K,L)

MU_FZF = zeros(2,L);

for k1 = 1:K
    if pilotIndex(k1)~=pilotIndex(k)
        for l = 1:L
            MU_FZF(1,l) = MU_FZF(1,l) + alpha^2*p_uplink(k1)*conj(A_FZF(l,k,1))*conj(A_FZF(l,k,1))*(beta_matrix(l,k1)-gamma_matrix(l,k1))/((N-tau_p)*theta_matrix(l,k));
            MU_FZF(2,l) = MU_FZF(2,l) + alpha^2*p_uplink(k1)*conj(A_FZF(l,k,2))*conj(A_FZF(l,k,2))*(beta_matrix(l,k1)-gamma_matrix(l,k1))/((N-tau_p)*theta_matrix(l,k));
        end
    end
end
end

function GN_FZF = computeGN_FZF(alpha, A_FZF, N, tau_p, theta_matrix, k,K,L)

GN_FZF = zeros(2,L);
for l = 1:L
    GN_FZF(1,l) = alpha^2*conj(A_FZF(l,k,1))*conj(A_FZF(l,k,1))/((N-tau_p)*theta_matrix(l,k));
    GN_FZF(2,l) = alpha^2*conj(A_FZF(l,k,2))*conj(A_FZF(l,k,2))/((N-tau_p)*theta_matrix(l,k));
end
end

function QN_FZF = computeQN_FZF(alpha, A_FZF, N, tau_p, theta_matrix, k, beta_matrix,K,L,p_uplink)

QN_FZF = zeros(2,L);
for l = 1:L
    QN_FZF_1 = 0;
    for k1 = 1:K
        QN_FZF_1 = QN_FZF_1 + p_uplink(k1)*beta_matrix(l,k1);
    end
    QN_FZF(1,l) =  alpha^3*(1 - alpha)*conj(A_FZF(l,k,1))*conj(A_FZF(l,k,1))*(QN_FZF_1+1)/((N-tau_p)*theta_matrix(l,k));
    QN_FZF(2,l) =  alpha^3*(1 - alpha)*conj(A_FZF(l,k,2))*conj(A_FZF(l,k,2))*(QN_FZF_1+1)/((N-tau_p)*theta_matrix(l,k));
end
end
