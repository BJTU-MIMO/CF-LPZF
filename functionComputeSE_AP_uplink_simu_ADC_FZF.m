function [DS, BU, PC, MU, GN, QN, A] = functionComputeSE_AP_uplink_simu_ADC_FZF(Hhat,H,R,C,nbrOfRealizations,N,K,L,D,p_uplink,alpha,pilotIndex,H_bar,tau_p)
%% Prepare to save simulation results
DS = zeros(K,2,3);
BU = zeros(K,2,3);
PC = zeros(K,2,3);
MU = zeros(K,2,3);
GN = zeros(K,2,3);
QN = zeros(K,2,3);

%%
G_ki  = zeros(L,K,K,3,nbrOfRealizations)  ;
G_kiki= zeros(L,L,K,K,3,nbrOfRealizations);

f2_k  = zeros(L,K,3,nbrOfRealizations)    ;
f2_kk = zeros(L,L,K,3,nbrOfRealizations)  ;

f1_k  = zeros(L,K,3,nbrOfRealizations)    ;
f1_kk = zeros(L,L,K,3,nbrOfRealizations)  ;

E_numer       =  zeros(L,K,3)      ;
E_denumer_gki =  zeros(L,L,K,K,3)  ;
E_denumer_f2k =  zeros(L,L,K,3)    ;
E_denumer_f1k =  zeros(L,L,K,3)    ;

A = zeros(L,K,3);

%combining vector
V_MR = zeros(N,L,K);
V_LPMMSE = zeros(N,L,K);
V_FZF = zeros(N,L,K);
%quantization error
R_q = zeros(N,N,L);

eyetau_p = eye(tau_p);
eyek = eye(K);

%% Go through all channel realizations


for k = 1:K
    
    for n = 1:nbrOfRealizations
        
        for l = 1:L
            
            %R_q Matrix
            for k1 = 1:K
                R_q(:,:,l) = R_q(:,:,l)+p_uplink(k1)*R(:,:,l,k1);
            end
            R_q(:,:,l) = alpha*(1-alpha)*diag(diag(squeeze(R_q(:,:,l))+eye(N)));
            
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
            V_FZF(:,l,k) = squeeze(H_bar((l-1)*N+1:l*N,n,:))/(squeeze(H_bar((l-1)*N+1:l*N,n,:))'*squeeze(H_bar((l-1)*N+1:l*N,n,:)))*eyetau_p(:,pilotIndex(k));
%             V_FZF(:,l,k) = squeeze(Hhat((l-1)*N+1:l*N,n,:))/(squeeze(Hhat((l-1)*N+1:l*N,n,:))'*squeeze(Hhat((l-1)*N+1:l*N,n,:)))*eyek(:,k);
%             V_FZF(:,l,k) = squeeze(H_bar((l-1)*N+1:l*N,n,k))/(squeeze(H_bar((l-1)*N+1:l*N,n,:))'*squeeze(H_bar((l-1)*N+1:l*N,n,:)));


            for i1 = 1:K
                %g_ki
                G_ki(l,k,i1,1,n) = squeeze(V_MR(:,l,k))'    *squeeze(D(:,:,l,k))*squeeze(H((l-1)*N+1:l*N,n,i1));
                G_ki(l,k,i1,2,n) = squeeze(V_LPMMSE(:,l,k))'*squeeze(D(:,:,l,k))*squeeze(H((l-1)*N+1:l*N,n,i1));
                G_ki(l,k,i1,3,n) = squeeze(V_FZF(:,l,k))'   *squeeze(H((l-1)*N+1:l*N,n,i1));
            end
            
            %f2_k
            f2_k(l,k,1,n)    = squeeze(V_MR(:,l,k))'    *squeeze(D(:,:,l,k))*chol(squeeze(R_q(:,:,l)))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1)) ;
            f2_k(l,k,2,n)    = squeeze(V_LPMMSE(:,l,k))'*squeeze(D(:,:,l,k))*chol(squeeze(R_q(:,:,l)))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1)) ;
            f2_k(l,k,3,n)    = squeeze(V_FZF(:,l,k))'   *chol(squeeze(R_q(:,:,l)))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1)) ;
            
            %f1_k
            f1_k(l,k,1,n)  = squeeze(V_MR(:,l,k))'      *squeeze(D(:,:,l,k))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1)) ;
            f1_k(l,k,2,n)  = squeeze(V_LPMMSE(:,l,k))'  *squeeze(D(:,:,l,k))*sqrt(0.5)*(randn(N,1) + 1i*randn(N,1)) ;
            f1_k(l,k,3,n)  = squeeze(V_FZF(:,l,k))'     *sqrt(0.5)*(randn(N,1) + 1i*randn(N,1)) ;
        end
        
        for i2 = 1:K
            G_kiki(:,:,k,i2,1,n) = squeeze(G_ki(:,k,i2,1,n))*squeeze(G_ki(:,k,i2,1,n))';
            G_kiki(:,:,k,i2,2,n) = squeeze(G_ki(:,k,i2,2,n))*squeeze(G_ki(:,k,i2,2,n))';
            G_kiki(:,:,k,i2,3,n) = squeeze(G_ki(:,k,i2,3,n))*squeeze(G_ki(:,k,i2,3,n))';
        end
        f2_kk(:,:,k,1,n) = squeeze(f2_k(:,k,1,n))*squeeze(f2_k(:,k,1,n))';
        f2_kk(:,:,k,2,n) = squeeze(f2_k(:,k,2,n))*squeeze(f2_k(:,k,2,n))';
        f2_kk(:,:,k,3,n) = squeeze(f2_k(:,k,3,n))*squeeze(f2_k(:,k,3,n))';
        f1_kk(:,:,k,1,n) = squeeze(f1_k(:,k,1,n))*squeeze(f1_k(:,k,1,n))';
        f1_kk(:,:,k,2,n) = squeeze(f1_k(:,k,2,n))*squeeze(f1_k(:,k,2,n))';
        f1_kk(:,:,k,3,n) = squeeze(f1_k(:,k,3,n))*squeeze(f1_k(:,k,3,n))';
    end
    
    %% expectation
    E_numer(:,k,1) = mean(G_ki(:,k,k,1,:),5);
    E_numer(:,k,2) = mean(G_ki(:,k,k,2,:),5);
    E_numer(:,k,3) = mean(G_ki(:,k,k,3,:),5);
    for i3 = 1:K
        E_denumer_gki(:,:,k,i3,1) = mean(G_kiki(:,:,k,i3,1,:),6);
        E_denumer_gki(:,:,k,i3,2) = mean(G_kiki(:,:,k,i3,2,:),6);
        E_denumer_gki(:,:,k,i3,3) = mean(G_kiki(:,:,k,i3,3,:),6);
    end
    E_denumer_f2k(:,:,k,1) = mean(f2_kk(:,:,k,1,:),5);
    E_denumer_f2k(:,:,k,2) = mean(f2_kk(:,:,k,2,:),5);
    E_denumer_f2k(:,:,k,3) = mean(f2_kk(:,:,k,3,:),5);
    E_denumer_f1k(:,:,k,1) = mean(f1_kk(:,:,k,1,:),5);
    E_denumer_f1k(:,:,k,2) = mean(f1_kk(:,:,k,2,:),5);
    E_denumer_f1k(:,:,k,3) = mean(f1_kk(:,:,k,3,:),5);
    
    %% Level 2
    a_k = (1/L) * ones(L,1);
    %MR
    DS(k,1,1) = abs(a_k'*alpha*sqrt(p_uplink(k))*E_numer(:,k,1))^2;
    BU(k,1,1) = a_k'*(alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,1) - alpha^2*p_uplink(k)*E_numer(:,k,1)*E_numer(:,k,1)')*a_k;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                PC(k,1,1) = PC(k,1,1) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,1))*a_k;
            end
        else
            MU(k,1,1) = MU(k,1,1) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,1))*a_k;
        end
    end
    GN(k,1,1) = alpha^2*a_k'*E_denumer_f1k(:,:,k,1)*a_k;
    QN(k,1,1) = a_k'*E_denumer_f2k(:,:,k,1)*a_k;
    
    %LPMMSE
    DS(k,1,2) = abs(a_k'*alpha*sqrt(p_uplink(k))*E_numer(:,k,2))^2;
    BU(k,1,2) = a_k'*(alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,2) - alpha^2*p_uplink(k)*E_numer(:,k,2)*E_numer(:,k,2)')*a_k;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
           if k1~=k
                PC(k,1,2) = PC(k,1,2) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,2))*a_k;
            end
        else
            MU(k,1,2) = MU(k,1,2) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,2))*a_k;
        end
    end
    GN(k,1,2) = alpha^2*a_k'*E_denumer_f1k(:,:,k,2)*a_k;
    QN(k,1,2) = a_k'*E_denumer_f2k(:,:,k,2)*a_k;
    
    %FZF
    DS(k,1,3) = abs(a_k'*alpha*sqrt(p_uplink(k))*E_numer(:,k,3))^2;
    BU(k,1,3) = a_k'*(alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,3) - alpha^2*p_uplink(k)*E_numer(:,k,3)*E_numer(:,k,3)')*a_k;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
           if k1~=k
                PC(k,1,3) = PC(k,1,3) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,3))*a_k;
            end
        else
            MU(k,1,3) = MU(k,1,3) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,3))*a_k;
        end
    end
    GN(k,1,3) = alpha^2*a_k'*E_denumer_f1k(:,:,k,3)*a_k;
    QN(k,1,3) = a_k'*E_denumer_f2k(:,:,k,3)*a_k;
    
   %% A
    A2 = alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,1);
    A3 = 0;
    A4 = 0;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                A3 = A3 + alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,1);
            end
        end
        if pilotIndex(k1)~=pilotIndex(k)
            A4 = A4 + alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,1);
        end
    end
    A5 = alpha^2*E_denumer_f1k(:,:,k,1);
    A6 = E_denumer_f2k(:,:,k,1);
    temp_1 = (A2+A3+A4+A5+A6-alpha^2*p_uplink(k)*E_numer(:,k,1)*squeeze(E_numer(:,k,1)'));
    temp_1(all(temp_1==0,2),:)=[];
%     temp_1(:,all(temp_1==0,1))=[];
    temp_2 = E_numer(:,k,1);
    temp_2(find(temp_2==0))=[];
    A(:,k,1) = sqrt(alpha^2*p_uplink(k))*((temp_1)\temp_2);
    
    A2 = alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,2);
    A3 = 0;
    A4 = 0;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                A3 = A3 + alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,2);
            end
        end
        if pilotIndex(k1)~=pilotIndex(k)
            A4 = A4 + alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,2);
        end
    end
    A5 = alpha^2*E_denumer_f1k(:,:,k,2);
    A6 = E_denumer_f2k(:,:,k,2);
    temp_1 = (A2+A3+A4+A5+A6-alpha^2*p_uplink(k)*E_numer(:,k,2)*squeeze(E_numer(:,k,2)'));
    temp_1(all(temp_1==0,2),:)=[];
%     temp_1(:,all(temp_1==0,1))=[];
    temp_2 = E_numer(:,k,2);
    temp_2(find(temp_2==0))=[];
    A(:,k,2) = sqrt(alpha^2*p_uplink(k))*((temp_1)\temp_2);
    
    A2 = alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,3);
    A3 = 0;
    A4 = 0;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                A3 = A3 + alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,3);
            end
        end
        if pilotIndex(k1)~=pilotIndex(k)
            A4 = A4 + alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,3);
        end
    end
    A5 = alpha^2*E_denumer_f1k(:,:,k,3);
    A6 = E_denumer_f2k(:,:,k,3);
    temp_1 = (A2+A3+A4+A5+A6-alpha^2*p_uplink(k)*E_numer(:,k,3)*squeeze(E_numer(:,k,3)'));
    temp_1(all(temp_1==0,2),:)=[];
%     temp_1(:,all(temp_1==0,1))=[];
    temp_2 = E_numer(:,k,3);
    temp_2(find(temp_2==0))=[];
    A(:,k,3) = sqrt(alpha^2*p_uplink(k))*((temp_1)\temp_2);
    
    %% Level 3
    a_k = squeeze(A(:,k,1));
    %MR
    DS(k,2,1) = abs(a_k'*alpha*sqrt(p_uplink(k))*E_numer(:,k,1))^2;
    BU(k,2,1) = a_k'*(alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,1) - alpha^2*p_uplink(k)*E_numer(:,k,1)*E_numer(:,k,1)')*a_k;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
            if k1~=k
                PC(k,2,1) = PC(k,2,1) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,1))*a_k;
            end
        else
            MU(k,2,1) = MU(k,2,1) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,1))*a_k;
        end
    end
    GN(k,2,1) = alpha^2*a_k'*E_denumer_f1k(:,:,k,1)*a_k;
    QN(k,2,1) = a_k'*E_denumer_f2k(:,:,k,1)*a_k;
    
    a_k = squeeze(A(:,k,2));
    %LPMMSE
    DS(k,2,2) = abs(a_k'*alpha*sqrt(p_uplink(k))*E_numer(:,k,2))^2;
    BU(k,2,2) = a_k'*(alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,2) - alpha^2*p_uplink(k)*E_numer(:,k,2)*E_numer(:,k,2)')*a_k;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
           if k1~=k
                PC(k,2,2) = PC(k,2,2) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,2))*a_k;
            end
        else
            MU(k,2,2) = MU(k,2,2) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,2))*a_k;
        end
    end
    GN(k,2,2) = alpha^2*a_k'*E_denumer_f1k(:,:,k,2)*a_k;
    QN(k,2,2) = a_k'*E_denumer_f2k(:,:,k,2)*a_k; 
    
    a_k = squeeze(A(:,k,3));
    %FZF
    DS(k,2,3) = abs(a_k'*alpha*sqrt(p_uplink(k))*E_numer(:,k,3))^2;
    BU(k,2,3) = a_k'*(alpha^2*p_uplink(k)*E_denumer_gki(:,:,k,k,3) - alpha^2*p_uplink(k)*E_numer(:,k,3)*E_numer(:,k,3)')*a_k;
    for k1 = 1:K
        if pilotIndex(k1)==pilotIndex(k)
           if k1~=k
                PC(k,2,3) = PC(k,2,3) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,3))*a_k;
            end
        else
            MU(k,2,3) = MU(k,2,3) + a_k'*(alpha^2*p_uplink(k1)*E_denumer_gki(:,:,k,k1,3))*a_k;
        end
    end
    GN(k,2,3) = alpha^2*a_k'*E_denumer_f1k(:,:,k,3)*a_k;
    QN(k,2,3) = a_k'*E_denumer_f2k(:,:,k,3)*a_k; 
end

