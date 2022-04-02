
%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 100;

%Number of channel realizations per setup
nbrOfRealizations = 100;

%Number of APs in the cell-free network
L = 10;

%Number of UEs in the network
K = 4;

%Number of antennas per AP
N = 4;

%Length of the coherence block
tau_c = 200;

%Compute number of pilots per coherence block
tau_p = 3;

%Uplink transmit power per UE for pilot(mW)
p_pilot = 100;


%%
%D matrix
D=zeros(N,N,L,K);
%Prepare to save results
SE_AP_MR_simu = zeros(K,2,nbrOfSetups);
SE_AP_LPMMSE_simu = zeros(K,2,nbrOfSetups);
SE_AP_FZF_simu = zeros(K,2,nbrOfSetups);
SE_AP_LPZF_simu = zeros(K,2,nbrOfSetups);
SE_AP_MR_anal = zeros(K,2,nbrOfSetups);
SE_AP_FZF_anal = zeros(K,2,nbrOfSetups);
SE_AP_LPZF_anal = zeros(K,2,nbrOfSetups);

%
DS_simu = zeros(K,2,4,nbrOfSetups);
BU_simu = zeros(K,2,4,nbrOfSetups);
PC_simu = zeros(K,2,4,nbrOfSetups);
MU_simu = zeros(K,2,4,nbrOfSetups);
QN_simu = zeros(K,2,4,nbrOfSetups);
GN_simu = zeros(K,2,4,nbrOfSetups);
A_simu  = zeros(L,K,4,nbrOfSetups);

DS_anal = zeros(K,2,3,nbrOfSetups);
BU_anal = zeros(K,2,3,nbrOfSetups);
PC_anal = zeros(K,2,3,nbrOfSetups);
MU_anal = zeros(K,2,3,nbrOfSetups);
QN_anal = zeros(K,2,3,nbrOfSetups);
GN_anal = zeros(K,2,3,nbrOfSetups);
A_anal  = zeros(L,K,2,3,nbrOfSetups);

%%
%Uplink transmit power per UE for uplink(mW)
p_uplink = 100*ones(K,1);

%%
prelogFactor = (1-tau_p/tau_c);

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndex,D_orig] = generateSetup_new_narrow(L,K,N,tau_p,1);
    beta_matrix = db2pow(gainOverNoisedB);
    
    %resoluation
    b = 12;
    alpha = finda(b);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,Psi,C,H_bar,gamma_matrix, theta_matrix, c_matrix] = functionChannelEstimatesADC_FZF(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p_pilot,alpha);


    for k = 1:K
        for l = 1:L
             D(:,:,l,k) = eye(N);
        end
    end
    
    [S,W,Z,M_matrix,E_S_temp] = functionUEgrouping_2(L,K,beta_matrix,pilotIndex,tau_p);
    
    %% simu
    [DS_simu(:,:,:,n), BU_simu(:,:,:,n), PC_simu(:,:,:,n), MU_simu(:,:,:,n), GN_simu(:,:,:,n), QN_simu(:,:,:,n), A_simu(:,:,:,n)] = functionComputeSE_AP_uplink_simu_ADC_LPZF(Hhat,H,R,C,nbrOfRealizations,N,K,L,D,p_uplink,alpha,pilotIndex,H_bar,tau_p,S,W,Z,M_matrix,E_S_temp,c_matrix);
   for k = 1:K
        SE_AP_MR_simu(k,1,n)     = prelogFactor*real(log2(1 + DS_simu(k,1,1,n)/(BU_simu(k,1,1,n)+PC_simu(k,1,1,n)+MU_simu(k,1,1,n)+GN_simu(k,1,1,n)+QN_simu(k,1,1,n))));   
        SE_AP_MR_simu(k,2,n)     = prelogFactor*real(log2(1 + DS_simu(k,2,1,n)/(BU_simu(k,2,1,n)+PC_simu(k,2,1,n)+MU_simu(k,2,1,n)+GN_simu(k,2,1,n)+QN_simu(k,2,1,n))));   
        SE_AP_LPMMSE_simu(k,1,n) = prelogFactor*real(log2(1 + DS_simu(k,1,2,n)/(BU_simu(k,1,2,n)+PC_simu(k,1,2,n)+MU_simu(k,1,2,n)+GN_simu(k,1,2,n)+QN_simu(k,1,2,n))));   
        SE_AP_LPMMSE_simu(k,2,n) = prelogFactor*real(log2(1 + DS_simu(k,2,2,n)/(BU_simu(k,2,2,n)+PC_simu(k,2,2,n)+MU_simu(k,2,2,n)+GN_simu(k,2,2,n)+QN_simu(k,2,2,n))));   
        SE_AP_FZF_simu(k,1,n)    = prelogFactor*real(log2(1 + DS_simu(k,1,3,n)/(BU_simu(k,1,3,n)+PC_simu(k,1,3,n)+MU_simu(k,1,3,n)+GN_simu(k,1,3,n)+QN_simu(k,1,3,n))));   
        SE_AP_FZF_simu(k,2,n)    = prelogFactor*real(log2(1 + DS_simu(k,2,3,n)/(BU_simu(k,2,3,n)+PC_simu(k,2,3,n)+MU_simu(k,2,3,n)+GN_simu(k,2,3,n)+QN_simu(k,2,3,n))));   
        SE_AP_LPZF_simu(k,1,n)   = prelogFactor*real(log2(1 + DS_simu(k,1,4,n)/(BU_simu(k,1,4,n)+PC_simu(k,1,4,n)+MU_simu(k,1,4,n)+GN_simu(k,1,4,n)+QN_simu(k,1,4,n))));   
        SE_AP_LPZF_simu(k,2,n)   = prelogFactor*real(log2(1 + DS_simu(k,2,4,n)/(BU_simu(k,2,4,n)+PC_simu(k,2,4,n)+MU_simu(k,2,4,n)+GN_simu(k,2,4,n)+QN_simu(k,2,4,n))));   
    end  
    
    %% anal
    [DS_anal(:,:,:,n), BU_anal(:,:,:,n), PC_anal(:,:,:,n), MU_anal(:,:,:,n), GN_anal(:,:,:,n), QN_anal(:,:,:,n), A_anal(:,:,:,:,n)] = functionComputeSE_AP_uplink_anal_ADC_LPZF_final(Psi,R,K,L,D,p_uplink,alpha,pilotIndex, p_pilot,tau_p,N,c_matrix,gamma_matrix, beta_matrix, theta_matrix, S,W,Z,M_matrix,E_S_temp);

    for k = 1:K
        SE_AP_MR_anal(k,1,n)      = prelogFactor*real(log2(1 + DS_anal(k,1,1,n)/(BU_anal(k,1,1,n)+PC_anal(k,1,1,n)+MU_anal(k,1,1,n)+GN_anal(k,1,1,n)+QN_anal(k,1,1,n))));
        SE_AP_MR_anal(k,2,n)      = prelogFactor*real(log2(1 + DS_anal(k,2,1,n)/(BU_anal(k,2,1,n)+PC_anal(k,2,1,n)+MU_anal(k,2,1,n)+GN_anal(k,2,1,n)+QN_anal(k,2,1,n))));
        SE_AP_FZF_anal(k,1,n)     = prelogFactor*real(log2(1 + DS_anal(k,1,2,n)/(BU_anal(k,1,2,n)+PC_anal(k,1,2,n)+MU_anal(k,1,2,n)+GN_anal(k,1,2,n)+QN_anal(k,1,2,n))));
        SE_AP_FZF_anal(k,2,n)     = prelogFactor*real(log2(1 + DS_anal(k,2,2,n)/(BU_anal(k,2,2,n)+PC_anal(k,2,2,n)+MU_anal(k,2,2,n)+GN_anal(k,2,2,n)+QN_anal(k,2,2,n))));
        SE_AP_LPZF_anal(k,1,n)    = prelogFactor*real(log2(1 + DS_anal(k,1,3,n)/(BU_anal(k,1,3,n)+PC_anal(k,1,3,n)+MU_anal(k,1,3,n)+GN_anal(k,1,3,n)+QN_anal(k,1,3,n))));
        SE_AP_LPZF_anal(k,2,n)    = prelogFactor*real(log2(1 + DS_anal(k,2,3,n)/(BU_anal(k,2,3,n)+PC_anal(k,2,3,n)+MU_anal(k,2,3,n)+GN_anal(k,2,3,n)+QN_anal(k,2,3,n))));
    end
    
   %% 
    %Remove large matrices at the end of analyzing this setup
    clear PsiInv C H Hhat;
end

%% Plot results
figure;
hold on; box on;
plot(sort(reshape(DS_simu(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(DS_anal(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
legend({'L3 (MR) simulation','L3 (MR) analytical'},'Interpreter','Latex','Location','NorthWest');

figure;
hold on; box on;
% plot(sort(reshape(BU_simu(:,1,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% plot(sort(reshape(BU_anal(:,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(BU_simu(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(BU_anal(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);

figure;
hold on; box on;
% plot(sort(reshape(PC_simu(:,1,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% plot(sort(reshape(PC_anal(:,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
plot(sort(reshape(PC_simu(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(PC_anal(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
% figure;
% hold on; box on;
% % plot(sort(reshape(MU_simu(:,1,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% % plot(sort(reshape(MU_anal(:,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% plot(sort(reshape(MU_simu(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
% plot(sort(reshape(MU_anal(:,2,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% figure;
% hold on; box on;
% % plot(sort(reshape(QN_simu(:,1,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% % plot(sort(reshape(QN_anal(:,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% plot(sort(reshape(QN_simu(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
% plot(sort(reshape(QN_anal(:,2,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% figure;
% hold on; box on;
% % plot(sort(reshape(GN_simu(:,1,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
% % plot(sort(reshape(GN_anal(:,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% plot(sort(reshape(GN_simu(:,2,1,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
% plot(sort(reshape(GN_anal(:,2,:) ,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% 

% figure;
% hold on; box on; 
% plot(sort(reshape(SE_AP_MR_simu(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m-','LineWidth',2);
% plot(sort(reshape(SE_AP_MR_anal(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m--','LineWidth',2);
% plot(sort(reshape(SE_AP_LPMMSE_simu(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
% plot(sort(reshape(SE_AP_FZF_simu(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
% plot(sort(reshape(SE_AP_FZF_anal(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k--','LineWidth',2);
% plot(sort(reshape(SE_AP_LPZF_simu(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'y-','LineWidth',2);
% plot(sort(reshape(SE_AP_LPZF_anal(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'y--','LineWidth',2);
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'L2 (MR) simulation ','L2 (MR) analytical ','L2 (LP-MMSE) ','L2 (FZF) simulation ','L2 (FZF) analytical ','L2 (LPZF) simulation ','L2 (LPZF) analytical '},'Interpreter','Latex','Location','NorthWest');

figure;
hold on; box on; 
plot(sort(reshape(SE_AP_MR_simu(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m-','LineWidth',2);
plot(sort(reshape(SE_AP_MR_anal(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m--','LineWidth',2);
plot(sort(reshape(SE_AP_LPMMSE_simu(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_AP_FZF_simu(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(reshape(SE_AP_FZF_anal(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k--','LineWidth',2);
plot(sort(reshape(SE_AP_LPZF_simu(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_AP_LPZF_anal(:,2,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r--','LineWidth',2);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'L3 (MR) simulation ','L3 (MR) analytical ','L3 (LP-MMSE) ','L3 (FZF) simulation ','L3 (FZF) analytical ','L3 (LPZF) simulation ','L3 (LPZF) analytical '},'Interpreter','Latex','Location','NorthWest');

