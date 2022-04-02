%% Define simulation setup

%Number of setups with random UE locations
nbrOfSetups = 200;

%Number of channel realizations per setup
nbrOfRealizations = 200;

%Number of APs in the cell-free network
L = 25;

%Number of UEs in the network
K = 10;

%Number of antennas per AP
N = 32;

%Length of the coherence block
tau_c = 200;

%Compute number of pilots per coherence block
tau_p = 7;

%Uplink transmit power per UE for pilot(mW)
p_pilot = 100;


%%
%D matrix
D=zeros(N,N,L,K);
%Prepare to save results
SE_AP_MR_tot = zeros(K,2,nbrOfSetups);
SE_AP_MMSE_tot = zeros(K,2,nbrOfSetups);
SE_AP_FZF_tot = zeros(K,2,nbrOfSetups);
SE_AP_RZF_tot = zeros(K,2,nbrOfSetups);
SE_AP_LPZF_tot = zeros(K,2,nbrOfSetups);
SE_AP_LPPZF_tot = zeros(K,2,nbrOfSetups);
%%
%Uplink transmit power per UE for uplink(mW)
p_uplink = 100*ones(K,1);
p_regu = p_uplink;
p_max = 100*K;

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    tic;
    %Generate one setup with UEs at random locations
   [gainOverNoisedB,R,pilotIndex,D_orig] = generateSetup_new(L,K,N,tau_p,1);
    beta_matrix = db2pow(gainOverNoisedB);
    
    %resoluation
    b = 12;
    alpha = finda(b);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,Psi,C,H_bar,gamma_matrix, theta_matrix, c_matrix] = functionChannelEstimatesADC_FZF(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p_pilot,alpha); 
    
    %D matrix
    for k = 1:K
        for l = 1:L
            D(:,:,l,k) = eye(N);
        end
    end
    
    [S,W,Z,M_matrix,E_S_temp] = functionUEgrouping(L,K,beta_matrix,pilotIndex,tau_p);
    
    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
   [SE_MR,SE_MMSE,SE_FZF,SE_RZF] = functionComputeSE_AP_uplink_lowADC_FZF_RZF_final(Hhat,H,R,C,tau_c,tau_p,nbrOfRealizations,N,K,L,D,p_uplink,alpha,pilotIndex,H_bar,p_regu,c_matrix);
   [SE_LPZF,SE_LPPZF] = functionComputeSE_AP_uplink_lowADC_LPZF_LPPZF_final(Hhat,H,R,tau_c,tau_p,nbrOfRealizations,N,K,L,p_uplink,alpha,pilotIndex,H_bar,S,W,Z,M_matrix,E_S_temp,c_matrix);
    
    %Save SE values
    SE_AP_MR_tot(:,:,n) = SE_MR;
    SE_AP_MMSE_tot(:,:,n) = SE_MMSE;
    SE_AP_FZF_tot(:,:,n) = SE_FZF;
    SE_AP_RZF_tot(:,:,n) = SE_RZF;
    SE_AP_LPZF_tot(:,:,n) = SE_LPZF;
    SE_AP_LPPZF_tot(:,:,n) = SE_LPPZF;
 
   %% 
    %Remove large matrices at the end of analyzing this setup
    clear PsiInv C H Hhat;
    toc;
end

%% Plot results
figure;
hold on; box on;
plot(sort(reshape(SE_AP_MR_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_AP_MMSE_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);
plot(sort(reshape(SE_AP_FZF_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
plot(sort(reshape(SE_AP_RZF_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'m-','LineWidth',2);
plot(sort(reshape(SE_AP_LPZF_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'y-','LineWidth',2);
plot(sort(reshape(SE_AP_LPPZF_tot(:,1,:),[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'g-','LineWidth',5);
xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'L3 (MR)','L3 (L-MMSE) ','L3 (FZF) ','L3 (RZF) ','L3 (LPZF) ','L3 (LPPZF) '},'Interpreter','Latex','Location','NorthWest');
