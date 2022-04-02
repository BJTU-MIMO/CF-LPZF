function [Hhat,H, H_bar,gamma_matrix, theta_matrix, c_matrix] = functionChannelEstimates_mRZF(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p,beta_matrix)
%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(L*N,nbrOfRealizations,K)+1i*randn(L*N,nbrOfRealizations,K));

B      = zeros(size(R))    ;
C      = zeros(size(R))    ;
Psi    = zeros(N,N,L,tau_p);

c_matrix = zeros(L,K);
gamma_matrix = zeros(L,K);
theta_matrix = zeros(L,K);

%Fi_matrix
Fi_matrix = sqrt(tau_p)*eye(tau_p);

%Go through all channels and apply the spatial correlation matrices
for l = 1:L
    
    for k = 1:K
        
        %Apply correlation to the uncorrelated channel realizations
        Rsqrt = sqrtm(R(:,:,l,k));
        H((l-1)*N+1:l*N,:,k) = sqrt(0.5)*Rsqrt*H((l-1)*N+1:l*N,:,k);
        
    end
    
end


%% Perform channel estimation

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,L,tau_p) + 1i*randn(N,nbrOfRealizations,L,tau_p));


%Prepare to store results
Hhat = zeros(L*N,nbrOfRealizations,K);

H_bar = zeros(L*N,nbrOfRealizations,tau_p);


%Go through all APs
for l = 1:L
    
    yq  = zeros(N,tau_p,nbrOfRealizations);

    for n = 1:nbrOfRealizations
        
        for  k = 1:K
            
            yq(:,:,n) = yq(:,:,n) + sqrt(p)*H((l-1)*N+1:l*N,n,k)*Fi_matrix(:,pilotIndex(k))';
            
        end
        
        H_bar((l-1)*N+1:l*N,n,:) = yq(:,:,n)*Fi_matrix;
        
    end
    
    %Go through all pilots
    for t = 1:tau_p
        
        %Compute processed pilot signal for all UEs that use pilot t
        yp = sqrt(p*tau_p)*sum(H((l-1)*N+1:l*N,:,t==pilotIndex),3) + Np(:,:,l,t);
        
        %Compute the matrix that is inverted in the MMSE estimator
        Psi(:,:,l,t) = p*tau_p*sum(R(:,:,l,t==pilotIndex),4) + eyeN;
        
        %Go through all UEs that use pilot t
        for k = find(t==pilotIndex)'
            
            %Compute the MMSE estimate
            RPsi = R(:,:,l,k) / squeeze(Psi(:,:,l,t));
            Hhat((l-1)*N+1:l*N,:,k) = sqrt(p*tau_p)*RPsi*yp;
            
            %Compute the spatial correlation matrix of the estimate
            B(:,:,l,k) = p*tau_p*RPsi*R(:,:,l,k);
            %Compute sum of all estimation error correlation matrices at every BS
            C(:,:,l,k) = (R(:,:,l,k)-B(:,:,l,k));
            
            gamma_matrix(l,k) = sum(diag(B(:,:,l,k)))/N;
            
            
        end
        
    end
    
    for k = 1:K
        
        c_matrix_temp = 0; 
        
        for i = 1:K
            
            if pilotIndex(i) == pilotIndex(k)
                
                c_matrix_temp = c_matrix_temp + p*tau_p*beta_matrix(l,i);
                
            end
            
        end
        
        c_matrix(l,k) = sqrt(p*tau_p)*beta_matrix(l,k)/(c_matrix_temp + 1);
        theta_matrix(l,k) = gamma_matrix(l,k)/(c_matrix(l,k)^2);
        
    end

    
end

end

