function  beta = srt_beta( d )
%SRT_BETA Summary of this function goes here
%   Detailed explanation goes here
[~,J_range,T_range]=size(d);
beta=zeros(J_range+1,J_range,T_range);
lambda=3*10^8/(1.9*10^9);%f=1.9Ghz
ht=[100;10*ones(J_range,1)];
hr=10*ones(1,J_range);
for t=(1:T_range)
    beta(:,:,t)=(lambda./(4*pi.*d(:,:,t))).^2.*(2*sin(2*pi*ht*hr./(lambda.*d(:,:,t)))).^2;
end
beta(isnan(beta))=0;
end

