%%
clear
clc
delta_t=60;%单位：秒
T_range=4*3600/delta_t;%an hour
J_range=24;
N=3;
P_i_max=100;P_j_max=10;
BW=20*10^6;%bandwidth=20MHz
%pos_i=zeros(2,1);pos_j=ones(2,t_range,J_range);
[pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t);%rand(J_range,T_range)+1i*rand(J_range,T_range);
pos_i=[zeros(1,T_range);pos_j];
d=srt_dis(pos_i,pos_j);
beta=srt_beta(d);
beta(d>30*10^3)=0;
gamma=0*ones(J_range,J_range);
%%

% %% 用户数J_range变化
% iter_N=20;
% j_step=3;
% T_range=8*3600/delta_t;%an hour
% srt_E=zeros(2,iter_N);ref_E=zeros(1,iter_N);
% srt_unsatisfied=zeros(1,iter_N);ref_unsatisfied=zeros(1,iter_N);
% J_ranges=(3+(0:iter_N-1)*j_step);
% for i=(1:iter_N)
%     %用户数J_range变化
%     J_range=J_ranges(i);
%     fprintf('J_range=%d\n',J_range)
%     [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t);%rand(J_range,T_range)+1i*rand(J_range,T_range);
%     pos_i=[zeros(1,T_range);pos_j];
%     d=srt_dis(pos_i,pos_j);
%     beta=srt_beta(d);
%     beta(d>30*10^3)=0;
%     gamma=0*ones(J_range,J_range);
%     [srt_E(1,i),srt_E(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_E(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E_means=[srt_E(1,:)./J_ranges;srt_E(2,:)./J_ranges];
% ref_E_means=ref_E./J_ranges;
% save('.\experiments\J_ranges_3_60','J_ranges','srt_E_means','ref_E_means')
% plot(J_ranges,srt_E_means(1,:),'r',J_ranges,srt_E_means(2,:),'r--o',J_ranges,ref_E_means,'g');

% %% delta_t 变化
% delta_ts=[240,180,120,60,30,15];
% srt_E=zeros(2,length(delta_ts));ref_E=zeros(1,length(delta_ts));
% srt_unsatisfied=zeros(1,length(delta_ts));ref_unsatisfied=zeros(1,length(delta_ts));
% for i=(1:length(delta_ts))
%     delta_t=delta_ts(i);
%     fprintf('delta_t=%d\n',delta_t)
%     T_range=6*3600/delta_t;
%     [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t);%rand(J_range,T_range)+1i*rand(J_range,T_range);
%     pos_i=[zeros(1,T_range);pos_j];
%     d=srt_dis(pos_i,pos_j);
%     beta=srt_beta(d);
%     beta(d>30*10^3)=0;
%     gamma=0*ones(J_range,J_range);
%     [srt_E(1,i),srt_E(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_E(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E_means=[srt_E(1,:)./J_range;srt_E(2,:)./J_range];
% ref_E_means=ref_E(1,:)./J_range;
% save('.\experiments\delta_ts_15_240','delta_ts','srt_E_means','ref_E_means')
% plot(delta_ts,srt_E_means(1,:),'r',delta_ts,srt_E_means(2,:),'r--o',delta_ts,ref_E_means(1,:),'g');

% %% T_range 变化
% divs=[1,2,3,4,5,6,12,30,60,120,240];
% srt_Es=zeros(2,length(divs));ref_Es=zeros(1,length(divs));
% srt_unsatisfied=zeros(1,length(divs));ref_unsatisfied=zeros(1,length(divs));
% [ref_E,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
% for i=(1:length(divs))
%     div=divs(i);
%     C_unsatisfied=C_qos;
%     fprintf('div=%d\n',div)
%     T_range_partial=T_range/div;
%     ii=1;
%     while ii<=div && sum(C_unsatisfied)
%         pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%         pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%         beta_partial=beta(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%         gamma=0*ones(J_range,J_range);
%         srt_Es_partial=zeros(2,1);
%         [srt_Es_partial(1),srt_Es_partial(2),C_unsatisfied]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied,N,gamma,delta_t,BW);
%         srt_Es(:,i)=srt_Es(:,i)+srt_Es_partial;
%         ii=ii+1;
%     end
% end
% srt_Es_means=[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
% ref_Es_means=ref_E*ones(1,length(divs))./J_range;
% save('.\experiments\T_ranges_1_240','divs','srt_Es_means','ref_Es_means')
% plot(divs,srt_Es_means(1,:),'r',divs,srt_Es_means(2,:),'r--o',divs,ref_Es_means,'g');

%% C_qos 变化
multis=[0.1,0.2,0.5,0.8,1,1.2,1.5,2,3,4,5,6,7,8,9,10];
srt_Es=zeros(2,length(multis));ref_Es=zeros(1,length(multis));
srt_unsatisfied=zeros(1,length(multis));ref_unsatisfied=zeros(1,length(multis));
for i=(1:length(multis))
    multi=multis(i);
    C_qos_new=C_qos*multi;
    fprintf('multi=%2.2f\n',multi)
    [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,gamma,delta_t,BW);
    srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
    [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,gamma,delta_t,BW);
    ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
end
srt_Es_means=[srt_Es(1,:)./J_range./(multis);srt_Es(2,:)./J_range./(multis)];
ref_Es_means=ref_Es./J_range./(multis);
save('.\experiments\C_qos_0.1_10.mat','multis','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
plot(multis(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r',multis(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'r--o',multis(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g');
