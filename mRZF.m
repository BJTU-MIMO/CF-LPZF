%% Define simulation setup
%Number of UEs in the network
K = 2;

%Number of setups with random UE locations
nbrOfSetups = 500;

%Number of channel realizations per setup
nbrOfRealizations = 200;

%Number of APs in the cell-free network
L = 64;
%L = 200;

%Number of antennas per AP
N = K;

%Length of the coherence block
tau_c = 200;

%Compute number of pilots per coherence block
tau_p = K/2;

%Uplink transmit power per UE for pilot(mW)
p_pilot = 100;


%%
%D matrix
D=zeros(N,N,L,K);
%Prepare to save results
SE_AP_mRZF_simu_tot = zeros(K,nbrOfSetups);
SE_AP_mRZF_anal_tot = zeros(K,nbrOfSetups);

%%
%Uplink transmit power per UE for uplink(mW)
p_uplink = 100*ones(K,1);
p_max = 100;

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    tic;
    %Generate one setup with UEs at random locations
   [gainOverNoisedB,R,pilotIndex,D_orig] = generateSetup_new(L,K,N,tau_p,1);
    beta_matrix = db2pow(gainOverNoisedB);
    
    %Generate channel realizations, channel estimates, and estimation
    %error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,H_bar,gamma_matrix,theta_matrix,c_matrix] = functionChannelEstimates_mRZF(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p_pilot,beta_matrix); 

    %Compute SE for the Cell-free mMIMO system with Monte Carlo simulations
   [SE_mRZF_simu, SE_mRZF_anal] = functionComputeSE_AP_uplink_mRZF_0106(Hhat,H,H_bar,tau_c,tau_p,nbrOfRealizations,N,K,L,p_uplink,pilotIndex,beta_matrix, gamma_matrix,theta_matrix,c_matrix);
   
    %Save SE values
    SE_AP_mRZF_simu_tot(:,n) = SE_mRZF_simu;
    SE_AP_mRZF_anal_tot(:,n) = SE_mRZF_anal;

 
   %% 
    %Remove large matrices at the end of analyzing this setup
    clear PsiInv C H Hhat;
    toc;
end

%% Plot results
figure;
hold on; box on;
plot(sort(reshape(SE_AP_mRZF_simu_tot,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'b-','LineWidth',2);
plot(sort(reshape(SE_AP_mRZF_anal_tot,[K*nbrOfSetups 1])),linspace(0,1,K*nbrOfSetups),'r-','LineWidth',2);

xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'simu','anal'},'Interpreter','Latex','Location','NorthWest');
