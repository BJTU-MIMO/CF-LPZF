function [ p_uplink ] = fractionalPowerControl(beta_matrix, p_max, K)
p_uplink = zeros(K,1); 
p_uplink_temp = 1./sum(beta_matrix);
for k = 1:K
    p_uplink(k,1) =  p_max*(p_uplink_temp(1,k)/max(p_uplink_temp));
end
end
% p_uplink = zeros(K,1); 
% for k = 1:K
%     p_uplink(k,1) =  K*p_max*(sum(beta_matrix(:,k))/sum(sum(beta_matrix)));
% end
% end

