function [SE_LPZF,SE_LPPWZF,SE_LPPSZF] = functionComputeSE_AP_uplink_lowADC_LPZF_LPPZF_final_PWandPS1123(Hhat,H,R,tau_c,tau_p,nbrOfRealizations,N,K,L,p_uplink,alpha,pilotIndex,H_bar,S,W,Z,M,E_S_temp,E_W_temp,c_matrix)
%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_LPZF = zeros(K,2);
SE_LPPWZF = zeros(K,2);
SE_LPPSZF = zeros(K,2);

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

signal_LPPWZF_level23 = zeros(L,K);
scaling_LPPWZF_level23 = zeros(L,K);
Gp_LPPWZF_level23 = zeros(L,L,K);
QN_LPPWZF_level23 = zeros(L,L,K);

signal_LPPSZF_level23 = zeros(L,K);
scaling_LPPSZF_level23 = zeros(L,K);
Gp_LPPSZF_level23 = zeros(L,L,K);
QN_LPPSZF_level23 = zeros(L,L,K);

V_LPZF = zeros(N,L,K);
V_LPPWZF = zeros(N,L,K);
V_LPPSZF = zeros(N,L,K);
%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Levels 2-3
    gp_LPZF_level23 = zeros(L,K,K);
    gp_LPPWZF_level23 = zeros(L,K,K);
    gp_LPPSZF_level23 = zeros(L,K,K);
    qn_LPZF_level23 = zeros(L,K);
    qn_LPPWZF_level23 = zeros(L,K);
    qn_LPPSZF_level23 = zeros(L,K);
    
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
            V_LPPWZF(:,l_Z_k,k) = c_matrix(Z_k_ind(l_Z_k),k)*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l/(E_S_l'*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))'*squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l)*epso;
            
%             E_W_l_ind = find(E_W_temp(Z_k_ind(l_Z_k),:)~=0);
%             if sum(E_W_l_ind)~=0
%                 tau_W_l = length(E_W_l_ind);
%                 E_W_l = zeros(tau_p,tau_W_l);
%                 
%                 for r = 1:tau_W_l
%                     E_W_l(:,r) = eyetau_p(:,E_W_l_ind(r));
%                 end
%                 B_matrix_temp = squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_W_l;
%                 B_matrix = eyeN - B_matrix_temp/(B_matrix_temp'*B_matrix_temp)*B_matrix_temp';
%             else
%                 B_matrix = eye(N);
%             end
%             
%             V_LPPSZF(:,l_Z_k,k) = c_matrix(Z_k_ind(l_Z_k),k)*B_matrix*(squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))*E_S_l/(E_S_l'*(squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:))')*(squeeze(H_bar((Z_k_ind(l_Z_k)-1)*N+1:Z_k_ind(l_Z_k)*N,n,:)))*E_S_l)*epso);
            V_LPPSZF(:,l_Z_k,k) = V_LPZF(:,l_Z_k,k);
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
            V_LPPWZF(:,Z_k_ind_length + l_M_k,k) = B_matrix*squeeze(Hhat((M_k_ind(l_M_k)-1)*N+1:M_k_ind(l_M_k)*N,n,k));
            V_LPPSZF(:,Z_k_ind_length + l_M_k,k) = V_LPZF(:,Z_k_ind_length + l_M_k,k);
            
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
            
            
            v_LPPWZF = V_LPPWZF(:,l,k);
            h_LPPWZF_i = reshape(h_group(:,l,k,:),[N K]);
            signal_LPPWZF_level23(l,k) = signal_LPPWZF_level23(l,k) + sqrt(p_uplink(k))*alpha*(v_LPPWZF'*h_LPPWZF_i(:,k))/nbrOfRealizations;
            gp_LPPWZF_level23(l,:,k) = gp_LPPWZF_level23(l,:,k) + alpha*(v_LPPWZF'*h_LPPWZF_i)*Dp12;
            if M(k,l) == 1
                for i_LPPWZF = 1:K
                    if Z(i_LPPWZF,l) == 1
                        gp_LPPWZF_level23(l,i_LPPWZF,k) = 0;
                    end
                end
            end
            scaling_LPPWZF_level23(l,k) = scaling_LPPWZF_level23(l,k) + alpha^2*norm(v_LPPWZF).^2/nbrOfRealizations;
            qn_LPPWZF_level23(l,k) = qn_LPPWZF_level23(l,k) + (v_LPPWZF'*N_q);

            v_LPPSZF = V_LPPSZF(:,l,k);
            h_LPPSZF_i = reshape(h_group(:,l,k,:),[N K]);
            signal_LPPSZF_level23(l,k) = signal_LPPSZF_level23(l,k) + sqrt(p_uplink(k))*alpha*(v_LPPSZF'*h_LPPSZF_i(:,k))/nbrOfRealizations;
            gp_LPPSZF_level23(l,:,k) = gp_LPPSZF_level23(l,:,k) + alpha*(v_LPPSZF'*h_LPPSZF_i)*Dp12;
            if Z(k,l) == 1
                for i_LPPSZF = 1:K
                    if M(i_LPPSZF,l) == 1
                        gp_LPPSZF_level23(l,i_LPPSZF,k) = 0;
                    end
                end
            end
            scaling_LPPSZF_level23(l,k) = scaling_LPPSZF_level23(l,k) + alpha^2*norm(v_LPPSZF).^2/nbrOfRealizations;
            qn_LPPSZF_level23(l,k) = qn_LPPSZF_level23(l,k) + (v_LPPSZF'*N_q);
            
        end
        
    end
    
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        Gp_LPZF_level23(:,:,k)  = Gp_LPZF_level23(:,:,k)  + gp_LPZF_level23(:,:,k) *gp_LPZF_level23(:,:,k)' /nbrOfRealizations;
        Gp_LPPWZF_level23(:,:,k) = Gp_LPPWZF_level23(:,:,k) + gp_LPPWZF_level23(:,:,k)*gp_LPPWZF_level23(:,:,k)'/nbrOfRealizations;
        Gp_LPPSZF_level23(:,:,k) = Gp_LPPSZF_level23(:,:,k) + gp_LPPSZF_level23(:,:,k)*gp_LPPSZF_level23(:,:,k)'/nbrOfRealizations;
  
        QN_LPZF_level23(:,:,k)  = qn_LPZF_level23(:,k) *qn_LPZF_level23(:,k)' /nbrOfRealizations;
        QN_LPPWZF_level23(:,:,k) = qn_LPPWZF_level23(:,k)*qn_LPPWZF_level23(:,k)'/nbrOfRealizations;
        QN_LPPSZF_level23(:,:,k) = qn_LPPSZF_level23(:,k)*qn_LPPSZF_level23(:,k)'/nbrOfRealizations;
        
    end
    
end

%Compute SE for Level 2 and Level 3
for k = 1:K
    
    b = signal_LPZF_level23(:,k);
    A = Gp_LPZF_level23(:,:,k) +  QN_LPZF_level23(:,:,k) + diag(scaling_LPZF_level23(:,k));
    A = A - (b*b');
    SE_LPZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));
    SE_LPZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));
    
    b = signal_LPPWZF_level23(:,k);
    A = Gp_LPPWZF_level23(:,:,k) +  QN_LPPWZF_level23(:,:,k) + diag(scaling_LPPWZF_level23(:,k)) ;
    A = A - (b*b');
    SE_LPPWZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));
    SE_LPPWZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));
    
    b = signal_LPPSZF_level23(:,k);
    A = Gp_LPPSZF_level23(:,:,k) +  QN_LPPSZF_level23(:,:,k) + diag(scaling_LPPSZF_level23(:,k)) ;
    A = A - (b*b');
    SE_LPPSZF(k,1) = prelogFactor*real(log2(1+b'*(A\b)));
    SE_LPPSZF(k,2) = prelogFactor*real(log2(1+abs(mean(b)).^2 / mean(mean(A))));
       
end




