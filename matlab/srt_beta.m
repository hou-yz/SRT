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
beta(d>30*10^3)=0;
% for t=(1:T_range)
%     d1=reshape(d(1,:,t),J_range,1);%到基站距离
%     in_cell_index=d1<30*10^3;%基站服务范围内用户
%     [js,~]=find(~in_cell_index);
%     for i=(1:length(js))
%         j=js(i);
%         beta(j+1,:,t)=zeros(1,J_range,1);
%         beta(:,j,t)=zeros(J_range+1,1,1);
%     end
% end
end

