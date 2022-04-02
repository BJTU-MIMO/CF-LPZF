function [SE_MR,SE_MMSE,SE_FZF,SE_RZF] = functionComputeSE_AP_uplink_lowADC_FZF_RZF_final(Hhat,H,R,C,tau_c,tau_p,nbrOfRealizations,N,K,L,D,p_uplink,alpha,pilotIndex,H_bar,p_regu,c_matrix)
%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MR = zeros(K,2);
SE_MMSE = zeros(K,2);
SE_FZF = zeros(K,2);
SE_RZF = zeros(K,2);

%quantization error
R_q = zeros(N,N,L);
for l = 1:L
    %R_q Matrix
    for k1 = 1:K
        R_q(:,:,l) = R_q(:,:,l) + p_uplink(k1)*R(:,:,l,k1);
    end
    R_q(:,:,l) = alpha*(1-alpha)*diag(diag(squeeze(R_q(:,:,l))+eye(N)));
end


%Diagonal matrix with transmit powers and its square root
Dp12 = diag(sqrt(p_uplink));

regu_matrix = diag(1./p_regu);

eyetau_p = eye(tau_p);
eyeN = eye(N);
eyeK = eye(K);


%Prepare to save simulation results
signal_MR_level23 = zeros(L,K);
scaling_MR_level23 = zeros(L,K);
Gp_MR_level23 = zeros(L,L,K);
QN_MR_level23 = zeros(L,L,K);

signal_MMSE_level23 = zeros(L,K);
scaling_MMSE_level23 = zeros(L,K);
G_MMSE_level23 = zeros(L,L,K);
QN_MMSE_level23 = zeros(L,L,K);

signal_FZF_level23 = zeros(L,K);
scaling_FZF_level23 = zeros(L,K);
Gq_FZF_level23 = zeros(L,L,K);
QN_FZF_level23 = zeros(L,L,K);

signal_RZF_level23 = zeros(L,K);
scaling_RZF_level23 = zeros(L,K);
Gq_RZF_level23 = zeros(L,L,K);
QN_RZF_level23 = zeros(L,L,K);

V_RZF = zeros(N,L,K);
V_FZF = zeros(N,L,K);
V_LPMMSE = zeros(N,L,K);
V_MR = zeros(N,L,K);

%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Levels 2-3
    gp_MR_level23 = zeros(L,K,K);
    gp_MMSE_level23 = zeros(L,K,K);
    gp_FZF_level23 = zeros(L,K,K);
    gp_RZF_level23 = zeros(L,K,K);
    qn_MR_level23 = zeros(L,K);
    qn_MMSE_level23 = zeros(L,K);
    qn_FZF_level23 = zeros(L,K);
    qn_RZF_level23 = zeros(L,K);
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        
        
        %
        N_q = chol(squeeze(R_q(:,:,l)))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1));
        
        
        for k = 1:K
            %combining vector
            %MR
            V_MR(:,l,k) = squeeze(D(:,:,l,k))*squeeze(Hhat((l-1)*N+1:l*N,n,k));
            %LPMMSE
            V_LPMMSE_temp = 0;
            for i = 1:K
                V_LPMMSE_temp = V_LPMMSE_temp + alpha^2*(sum(diag(squeeze(D(:,:,l,i))))/N)*p_uplink(i)*(squeeze(Hhat((l-1)*N+1:l*N,n,i))*squeeze(Hhat((l-1)*N+1:l*N,n,i))'+ C(:,:,l,i));
            end
            V_LPMMSE(:,l,k) = squeeze(V_LPMMSE_temp + alpha^2*eye(N) + R_q(:,:,l))\squeeze(Hhat((l-1)*N+1:l*N,n,k))*alpha*p_uplink(k);
            %FZF
            V_FZF(:,l,k) = c_matrix(l,k)*squeeze(H_bar((l-1)*N+1:l*N,n,:))/(squeeze(H_bar((l-1)*N+1:l*N,n,:))'*squeeze(H_bar((l-1)*N+1:l*N,n,:)))*eyetau_p(:,pilotIndex(k));
            %RZF
            V_RZF(:,l,k) = squeeze(Hhat((l-1)*N+1:l*N,n,:))/(squeeze(Hhat((l-1)*N+1:l*N,n,:))'*squeeze(Hhat((l-1)*N+1:l*N,n,:)) + regu_matrix)*eyeK(:,k);
            V_RZF_test = squeeze(H_bar((l-1)*N+1:l*N,n,:))/(squeeze(H_bar((l-1)*N+1:l*N,n,:))'*squeeze(H_bar((l-1)*N+1:l*N,n,:)) + pinv(diag(p_uplink(1:tau_p))))*eyetau_p(:,pilotIndex(k));
           
            
        end
        
        %Go through all UEs
        for k = 1:K
            
            
            %%MR combining
            v = V_MR(:,l,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_MR_level23(l,k) = signal_MR_level23(l,k) + sqrt(p_uplink(k))*alpha*(v'*D(:,:,l,k)*Hallj(:,k))/nbrOfRealizations;
            gp_MR_level23(l,:,k) = gp_MR_level23(l,:,k) + alpha*(v'*D(:,:,l,k)*Hallj)*Dp12;
            scaling_MR_level23(l,k) = scaling_MR_level23(l,k) + alpha^2*norm(v).^2/nbrOfRealizations;
            qn_MR_level23(l,k) = qn_MR_level23(l,k) + (v'*D(:,:,l,k)*N_q);
            
            
            %%MMSE combining
            v = V_LPMMSE(:,l,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_MMSE_level23(l,k) = signal_MMSE_level23(l,k) + sqrt(p_uplink(k))*alpha*(v'*D(:,:,l,k)*Hallj(:,k))/nbrOfRealizations;
            gp_MMSE_level23(l,:,k) = gp_MMSE_level23(l,:,k) + alpha*(v'*D(:,:,l,k)*Hallj)*Dp12;
            scaling_MMSE_level23(l,k) = scaling_MMSE_level23(l,k) + alpha^2*norm(v).^2/nbrOfRealizations;
            qn_MMSE_level23(l,k) = qn_MMSE_level23(l,k) + (v'*D(:,:,l,k)*N_q);
            
            
            %%FZF combining
            v = V_FZF(:,l,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_FZF_level23(l,k) = signal_FZF_level23(l,k) + sqrt(p_uplink(k))*alpha*(v'*D(:,:,l,k)*Hallj(:,k))/nbrOfRealizations;
            gp_FZF_level23(l,:,k) = gp_FZF_level23(l,:,k) + alpha*(v'*D(:,:,l,k)*Hallj)*Dp12;
            scaling_FZF_level23(l,k) = scaling_FZF_level23(l,k) + alpha^2*norm(v).^2/nbrOfRealizations;
            qn_FZF_level23(l,k) = qn_FZF_level23(l,k) + (v'*D(:,:,l,k)*N_q);
            
            %%RZF combining
            v = V_RZF(:,l,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_RZF_level23(l,k) = signal_RZF_level23(l,k) + sqrt(p_uplink(k))*alpha*(v'*D(:,:,l,k)*Hallj(:,k))/nbrOfRealizations;
            gp_RZF_level23(l,:,k) = gp_RZF_level23(l,:,k) + alpha*(v'*D(:,:,l,k)*Hallj)*Dp12;
            scaling_RZF_level23(l,k) = scaling_RZF_level23(l,k) + alpha^2*norm(v).^2/nbrOfRealizations;
            qn_RZF_level23(l,k) = qn_RZF_level23(l,k) + (v'*D(:,:,l,k)*N_q);


        end
        
    end
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        Gp_MR_level23(:,:,k) = Gp_MR_level23(:,:,k) + gp_MR_level23(:,:,k)*gp_MR_level23(:,:,k)'/nbrOfRealizations;
        G_MMSE_level23(:,:,k) = G_MMSE_level23(:,:,k) + gp_MMSE_level23(:,:,k)*gp_MMSE_level23(:,:,k)'/nbrOfRealizations;
        Gq_FZF_level23(:,:,k) = Gq_FZF_level23(:,:,k) + gp_FZF_level23(:,:,k)*gp_FZF_level23(:,:,k)'/nbrOfRealizations;
        Gq_RZF_level23(:,:,k) = Gq_RZF_level23(:,:,k) + gp_RZF_level23(:,:,k)*gp_RZF_level23(:,:,k)'/nbrOfRealizations;
        
        QN_MR_level23(:,:,k) = qn_MR_level23(:,k)*qn_MR_level23(:,k)'/nbrOfRealizations;
        QN_MMSE_level23(:,:,k) = qn_MMSE_level23(:,k)*qn_MMSE_level23(:,k)'/nbrOfRealizations;
        QN_FZF_level23(:,:,k) = qn_FZF_level23(:,k)*qn_FZF_level23(:,k)'/nbrOfRealizations;
        QN_RZF_level23(:,:,k) = qn_RZF_level23(:,k)*qn_RZF_level23(:,k)'/nbrOfRealizations;
        
    end
    
end

%Compute SE for Level 2 and Level 3
for k = 1:K
    %With MR combining
    b = signal_MR_level23(:,k);
    A = Gp_MR_level23(:,:,k) +  QN_MR_level23(:,:,k) + diag(scaling_MR_level23(:,k));
    A = A - (b*b');
    SE_MR(k,1) = prelogFactor*real(log2(1+b'*(A\b)));   
    SE_MR(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));   

    %With L-MMSE combining
    b = signal_MMSE_level23(:,k);
    A = G_MMSE_level23(:,:,k) +  QN_MMSE_level23(:,:,k) + diag(scaling_MMSE_level23(:,k)) ;
    A = A - (b*b');
    SE_MMSE(k,1) = prelogFactor*real(log2(1+b'*(A\b)));  
    SE_MMSE(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));   
    
   %With FZF combining
    b = signal_FZF_level23(:,k);
    A = Gq_FZF_level23(:,:,k) +  QN_FZF_level23(:,:,k) + diag(scaling_FZF_level23(:,k)) ;
    A = A - (b*b');
    SE_FZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));  
    SE_FZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));  
    
   %With RZF combining
    b = signal_RZF_level23(:,k);
    A = Gq_RZF_level23(:,:,k) +  QN_RZF_level23(:,:,k) + diag(scaling_RZF_level23(:,k)) ;
    A = A - (b*b');
    SE_RZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));  
    SE_RZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));   

end




