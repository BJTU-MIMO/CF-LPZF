function [S,W,Z,M,E_S_temp,E_W_temp] = functionUEgrouping_1123(L,K,beta_matrix,pilotIndex,tau_p)
S = zeros(L,K);
W = zeros(L,K);

Z = zeros(K,L);
M = zeros(K,L);

E_S_temp = zeros(L,tau_p);
E_W_temp = zeros(L,tau_p);

beta_sele_matrix = zeros(L,K);

thres = 0.85;

for l = 1:L
    [beta_des, beta_des_ind] = sort(beta_matrix(l,:),'descend');
    temp1 = 0;
    for k = 1:K
        temp1 = temp1 + beta_des(k)/sum(beta_matrix(l,:));
        if temp1 > thres && k~=1
            break;
        else
            beta_sele_matrix(l,beta_des_ind(k)) = 1;
        end
    end
end

for l = 1:L
    for  k = 1:K
        if beta_sele_matrix(l,k)==1
            S(l,k) = 1;
            Z(k,l) = 1;
            for k1 = 1:K
                if pilotIndex(k1)==pilotIndex(k)
                    S(l,k1) = 1;
                    Z(k1,l) = 1;
                end
            end
            E_S_temp(l,pilotIndex(k)) = 1;
%         else
%             W(l,k) = 1;
%             M(k,l) = 1;
%             for k1 = 1:K
%                 if pilotIndex(k1)==pilotIndex(k)
%                     W(l,k1) = 1;
%                     M(k1,l) = 1;
%                 end
%             end
%             E_W_temp(l,pilotIndex(k)) = 1;
        end
    end
end

I_temp1 = ones(L,K);
W = I_temp1 - S;
I_temp2 = ones(K,L);
M = I_temp2 - Z;
I_temp3 = ones(L,tau_p);
E_W_temp = I_temp3 - E_S_temp;

end

