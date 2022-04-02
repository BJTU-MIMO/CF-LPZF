function [SE_LPZF,SE_LPPZF] = functionComputeSE_AP_uplink_lowADC_LPZF_LPPZF_final(Hhat,H,R,tau_c,tau_p,nbrOfRealizations,N,K,L,p_uplink,alpha,pilotIndex,H_bar,S,W,Z,M,E_S_temp,c_matrix)
%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_LPZF = zeros(K,2);
SE_LPPZF = zeros(K,2);

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

eyeN = eye(N);
eyetau_p = eye(tau_p);


%Prepare to save simulation results
signal_LPZF_level23 = zeros(L,K);
scaling_LPZF_level23 = zeros(L,K);
Gp_LPZF_level23 = zeros(L,L,K);
QN_LPZF_level23 = zeros(L,L,K);

signal_LPPZF_level23 = zeros(L,K);
scaling_LPPZF_level23 = zeros(L,K);
Gp_LPPZF_level23 = zeros(L,L,K);
QN_LPPZF_level23 = zeros(L,L,K);

V_LPZF = zeros(N,L,K);
V_LPPZF = zeros(N,L,K);
%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Levels 2-3
    gp_LPZF_level23 = zeros(L,K,K);
    gp_LPPZF_level23 = zeros(L,K,K);
    qn_LPZF_level23 = zeros(L,K);
    qn_LPPZF_level23 = zeros(L,K);
    
    h_group = zeros(N,L,K,K);
    
    for k = 1:K
        %LPZF
        [~,Z_k_ind] = find(Z(k,:) == 1);
        Z_k_ind_length = length(Z_k_ind);
        
        %strong AP
        for l_Z_k = 1:Z_k_ind_length
            
            E_S_l_ind = find(E_S_temp(Z_k_ind(l_Z_k),:)~=0);
            tau_S_l = length(E_S_l_ind);
            E_S_l = zeros(tau_p,tau_S_l);
            eyetau_S_l = eye(tau_S_l );
            
            for r = 1:tau_S_l
                E_S_l(:,r) = eyetau_p(:,E_S_l_ind(r));
            end
           epso = eyetau_S_l(:,find(E_S_l_ind == pilotIndex(k)));
            V_LPZF(:,l_Z_k,k)  = c_matrix(Z_k_ind(l_Z_k),k)*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l/(E_S_l'*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))'*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l)*epso;
            V_LPPZF(:,l_Z_k,k) = c_matrix(Z_k_ind(l_Z_k),k)*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l/(E_S_l'*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))'*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l)*epso;
        end
        
        %Weak AP
        [~,M_k_ind] = find(M(k,:) == 1);
        M_k_ind_length = length(M_k_ind);
        for  l_M_k = 1:M_k_ind_length
            
            E_S_l_ind = find(E_S_temp(M_k_ind(l_M_k),:)~=0);
            if sum(E_S_l_ind)~=0
                tau_S_l = length(E_S_l_ind);
                E_S_l = zeros(tau_p,tau_S_l);
                
                for r = 1:tau_S_l
                    E_S_l(:,r) = eyetau_p(:,E_S_l_ind(r));
                end
                B_matrix_temp = squeeze(H_bar((M_k_ind(l_M_k)-1)*N+1:M_k_ind(l_M_k)*N,n,:))*E_S_l;
                B_matrix = eyeN - B_matrix_temp/(B_matrix_temp'*B_matrix_temp)*B_matrix_temp';
            else
                B_matrix = eye(N);
            end
            
            V_LPZF(:,Z_k_ind_length + l_M_k,k)  = squeeze(Hhat((M_k_ind(l_M_k)-1)*N+1:M_k_ind(l_M_k)*N,n,k));
            V_LPPZF(:,Z_k_ind_length + l_M_k,k) = B_matrix*squeeze(Hhat((M_k_ind(l_M_k)-1)*N+1:M_k_ind(l_M_k)*N,n,k));
        end
        
    end    
    
    for k = 1:K
        
        %Level 2 and Level 3
        
        %LPZF
        ind_Z_l = find(Z(k,:)~=0);
        ind_M_l = find(M(k,:)~=0);
        
        for count1 = 1:length(ind_Z_l)
            for k_i = 1:K
                Halli = reshape(H(:,n,k_i),[N L]);
                h_group(:,count1,k,k_i) = Halli(:,ind_Z_l(count1));
            end
        end
        for count2 = 1:length(ind_M_l)
            for k_i = 1:K
                Halli = reshape(H(:,n,k_i),[N L]);
                h_group(:,(length(ind_Z_l)+count2),k,k_i) = Halli(:,ind_M_l(count2));
            end
        end
        
    end
    
    for k = 1:K
        
        for l = 1:L
            
            N_q = chol(squeeze(R_q(:,:,l)))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1));
            
            v_LPZF = V_LPZF(:,l,k);
            h_LPZF_i = reshape(h_group(:,l,k,:),[N K]);
            signal_LPZF_level23(l,k) = signal_LPZF_level23(l,k) + sqrt(p_uplink(k))*alpha*(v_LPZF'*h_LPZF_i(:,k))/nbrOfRealizations;
            gp_LPZF_level23(l,:,k) = gp_LPZF_level23(l,:,k) + alpha*(v_LPZF'*h_LPZF_i)*Dp12;
            scaling_LPZF_level23(l,k) = scaling_LPZF_level23(l,k) + alpha^2*norm(v_LPZF).^2/nbrOfRealizations;
            qn_LPZF_level23(l,k) = qn_LPZF_level23(l,k) + (v_LPZF'*N_q);
            
            
            v_LPPZF = V_LPPZF(:,l,k);
            h_LPPZF_i = reshape(h_group(:,l,k,:),[N K]);
            signal_LPPZF_level23(l,k) = signal_LPPZF_level23(l,k) + sqrt(p_uplink(k))*alpha*(v_LPPZF'*h_LPPZF_i(:,k))/nbrOfRealizations;
            gp_LPPZF_level23(l,:,k) = gp_LPPZF_level23(l,:,k) + alpha*(v_LPPZF'*h_LPPZF_i)*Dp12;
            if M(k,l) == 1
                for i_LPPZF = 1:K
                    if Z(i_LPPZF,l) == 1
                        gp_LPPZF_level23(l,i_LPPZF,k) = 0;
                    end
                end
            end
            scaling_LPPZF_level23(l,k) = scaling_LPPZF_level23(l,k) + alpha^2*norm(v_LPPZF).^2/nbrOfRealizations;
            qn_LPPZF_level23(l,k) = qn_LPPZF_level23(l,k) + (v_LPPZF'*N_q);
        end
        
    end
    
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        Gp_LPZF_level23(:,:,k)  = Gp_LPZF_level23(:,:,k)  + gp_LPZF_level23(:,:,k) *gp_LPZF_level23(:,:,k)' /nbrOfRealizations;
        Gp_LPPZF_level23(:,:,k) = Gp_LPPZF_level23(:,:,k) + gp_LPPZF_level23(:,:,k)*gp_LPPZF_level23(:,:,k)'/nbrOfRealizations;
        
        QN_LPZF_level23(:,:,k)  = qn_LPZF_level23(:,k) *qn_LPZF_level23(:,k)' /nbrOfRealizations;
        QN_LPPZF_level23(:,:,k) = qn_LPPZF_level23(:,k)*qn_LPPZF_level23(:,k)'/nbrOfRealizations;
        
    end
    
end

%Compute SE for Level 2 and Level 3
for k = 1:K
    
    b = signal_LPZF_level23(:,k);
    A = Gp_LPZF_level23(:,:,k) +  QN_LPZF_level23(:,:,k) + diag(scaling_LPZF_level23(:,k));
    A = A - (b*b');
    SE_LPZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));
    SE_LPZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));
    
    b = signal_LPPZF_level23(:,k);
    A = Gp_LPPZF_level23(:,:,k) +  QN_LPPZF_level23(:,:,k) + diag(scaling_LPPZF_level23(:,k)) ;
    A = A - (b*b');
    SE_LPPZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));
    SE_LPPZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));
    
end




