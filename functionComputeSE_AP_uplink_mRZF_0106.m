function [SE_mRZF_simu, SE_mRZF_anal] = functionComputeSE_AP_uplink_mRZF_0106(Hhat,H,H_bar,tau_c,tau_p,nbrOfRealizations,N,K,L,p_uplink,pilotIndex,beta_matrix,gamma_matrix,theta_matrix,c_matrix)
%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_mRZF_simu = zeros(K,1);
SE_mRZF_anal = zeros(K,1);

%% simulation
Dp = diag(sqrt(p_uplink));
eyeN = eye(N);

%Prepare to save simulation results
signal_mRZF = zeros(L,K);
scaling_mRZF = zeros(L,K);
Gp_mRZF = zeros(L,L,K);

V_mRZF = zeros(N,L,K);
alpha = 0.8;

for n = 1:nbrOfRealizations
    
    %Levels 2-3
    gp_mRZF = zeros(L,K,K);
    
    %Go through all APs
    for l = 1:L
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(l-1)*N:l*N,n,:),[N K]);
        
        for k = 1:K
            
            %combining vector
            %mRZF
            V_mRZF(:,l,k) = c_matrix(l,k)*((squeeze(H_bar((l-1)*N+1:l*N,n,:))*squeeze(H_bar((l-1)*N+1:l*N,n,:))' + N*alpha*eyeN)\squeeze(H_bar((l-1)*N+1:l*N,n,pilotIndex(k))));
        end
        
        %Go through all UEs
        for k = 1:K
            
            %mRZF combining
            v = V_mRZF(:,l,k); %Extract combining vector
            
            %Level 2 and Level 3
            signal_mRZF(l,k) = signal_mRZF(l,k) + sqrt(p_uplink(k))*(v'*Hallj(:,k))/nbrOfRealizations;
            gp_mRZF(l,:,k) = gp_mRZF(l,:,k) + (v'*Hallj)*Dp;
            scaling_mRZF(l,k) = scaling_mRZF(l,k) + norm(v).^2/nbrOfRealizations;
            
        end
        
    end
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        Gp_mRZF(:,:,k) = Gp_mRZF(:,:,k) + gp_mRZF(:,:,k)*gp_mRZF(:,:,k)'/nbrOfRealizations;
        
    end
    
end

%Compute SE for Level 2 and Level 3
for k = 1:K
    
    %With RZF combining
    b = signal_mRZF(:,k);
    A = Gp_mRZF(:,:,k) + diag(scaling_mRZF(:,k)) ;
    A = A - (b*b');
    SE_mRZF_simu(k,1) = prelogFactor*real(log2(1+abs(mean(b)).^2 / abs(mean(mean(A)))));
    
end

%% analytical
eyeN = eye(N);
eyetaup = eye(tau_p);

a = 1;
DS = zeros(K,1);
PC_temp = zeros(K,K);
PC = zeros(K,1);
UI = zeros(K,1);
GN = zeros(K,1);

Big_theta_matrix = zeros(N,N,tau_p,L);
e_matrix = (1/alpha)*ones(tau_p,L);
e_hat_matrix = zeros(tau_p,L);
Big_e_hat_matrix = zeros(tau_p,tau_p,L);
J_matrix = zeros(tau_p,tau_p,L);
w_matrix = zeros(tau_p,L);
Big_w_matrix = zeros(tau_p, tau_p,L);
T_matrix = zeros(N,N,L);
zeta_matrix = zeros(tau_p,L);

count_TE = 10;

for l = 1:L
    
    for tau_k = 1:tau_p
        
        Big_theta_matrix_temp = 0;
        
        for t = 1:K
            
            if pilotIndex(t) == tau_k
                
                 Big_theta_matrix_temp = Big_theta_matrix_temp + p_uplink(t)*tau_p*beta_matrix(l,t) ;
                
            end
            
        end
        
        Big_theta_matrix(:,:,tau_k,l) = (Big_theta_matrix_temp + 1)*eyeN;
        
    end
    
end

for l = 1:L
    
    for count = 1:1:count_TE
        
        T_l_temp = 0;
        
        for tau_T = 1:tau_p
            
            T_l_temp = T_l_temp + squeeze(Big_theta_matrix(:,:,tau_T,l))/(1+e_matrix(tau_T,l));
            
        end
        
        T_matrix(:,:,l) = pinv((1/N)*T_l_temp+alpha*eyeN);
        
        for tau_e = 1:tau_p
            
            e_matrix(tau_e,l) = (1/N)*trace(squeeze(Big_theta_matrix(:,:,tau_e,l))*squeeze(T_matrix(:,:,l)));
            
        end
        
    end
    
end

for l = 1:L
      
    for tau_i = 1:tau_p
        
        for tau_j = 1:tau_p
            
            J_matrix(tau_i,tau_j,l) = (1/N)*trace(squeeze(Big_theta_matrix(:,:,tau_i,l))*squeeze(T_matrix(:,:,l))*squeeze(Big_theta_matrix(:,:,tau_j,l))*squeeze(T_matrix(:,:,l)))/(N*(1+e_matrix(tau_j,l))^2);
            
            Big_w_matrix(tau_i,tau_j,l) = (1/N)*trace(squeeze(Big_theta_matrix(:,:,tau_i,l))*squeeze(T_matrix(:,:,l))*squeeze(Big_theta_matrix(:,:,tau_j,l))*squeeze(T_matrix(:,:,l)));
            
        end
        
        w_matrix (tau_i,l) = (1/N)*trace(squeeze(Big_theta_matrix(:,:,tau_i,l))*squeeze(T_matrix(:,:,l))*squeeze(T_matrix(:,:,l)));
        
    end
    
    e_hat_matrix(:,l) = squeeze((eyetaup - J_matrix(:,:,l))\w_matrix(:,l));
    
    for tau_k1 = 1:tau_p
        
        Big_e_hat_matrix(:,tau_k1,l) =  (eyetaup - J_matrix(:,:,l))\squeeze(Big_w_matrix(:,tau_k1,l));
        
    end
    
    for tau_k2 = 1:tau_p
        
        zeta_matrix(tau_k2,l) = (1/N)*e_hat_matrix(tau_k2,l)/(1+e_matrix(tau_k2,l))^2;
        
    end
    
    
end

for k = 1:K
    
    for l = 1:L
        
        DS(k) = DS(k) + a*c_matrix(l,k)^2*e_matrix(pilotIndex(k),l)/(1+e_matrix(pilotIndex(k),l));
        
    end
    
    DS(k) = p_uplink(k)*DS(k)^2;
    
    for i = 1:K
            
        if pilotIndex(i) == pilotIndex(k) && i~=k
            
            for l = 1:L
                
                PC_temp(k,i) = PC_temp(k,i) + a*c_matrix(l,k)*c_matrix(l,i)*e_matrix(pilotIndex(k),l)/(1+e_matrix(pilotIndex(k),l));
                
            end
            
        end  
        
        if pilotIndex(i) ~= pilotIndex(k)
            
            for l  = 1:L
                
                UI(k) = UI(k) + (1/N)*p_uplink(i)*a^2*c_matrix(l,k)^2*(c_matrix(l,i)^2*theta_matrix(l,i)*Big_e_hat_matrix(pilotIndex(i),pilotIndex(k),l) + (beta_matrix(l,i)-gamma_matrix(l,i))*e_hat_matrix(pilotIndex(i),l)*(1+e_matrix(pilotIndex(i),l))^2)/theta_matrix(l,i)/(1+e_matrix(pilotIndex(k),l))^2/(1+e_matrix(pilotIndex(i),l))^2;
                
            end
            
        end
        
    end
    
    PC(k) = p_uplink(k)*sum(PC_temp(k,:))^2;
    
    for l = 1:L
        
        GN(k) = GN(k) + a^2*c_matrix(l,k)^2*zeta_matrix(pilotIndex(k),l);
        
    end
    
    SE_mRZF_anal(k) = prelogFactor*abs(log2(1+DS(k)/(PC(k)+UI(k)+GN(k))));
    
end

end



