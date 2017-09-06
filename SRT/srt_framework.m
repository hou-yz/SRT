%%
clear
clc
%航线1，航线2，航线3 发船间隔
interval=15*60;%单位：s
%航线1，航线2，航线3 占比
ratio=[0,1/2,1/2];
delta_t=60;%单位：秒
T_range=6*3600/delta_t;%an hour
J_range=(T_range*delta_t-2*3600)/interval*2;
N=1;
P_i_max=100;P_j_max=10;
BW=2*10^6;%bandwidth=2MHz

[pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
pos_i=[zeros(1,T_range);pos_j];
d=srt_dis(pos_i,pos_j);
beta=srt_beta(d);
gamma=0*ones(J_range,J_range);

% [srt_E_1,srt_E_2,C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
% srt_unsatisfied=(sum(C_unsatisfied)>0);
% [ref_Es,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
% ref_unsatisfied=(sum(C_unsatisfied)>0);

%% delta_t 变化
% delta_ts=[240,180,120,60,30,15];
% srt_Es=zeros(2,length(delta_ts));ref_Es=zeros(1,length(delta_ts));
% srt_unsatisfied=zeros(1,length(delta_ts));ref_unsatisfied=zeros(1,length(delta_ts));
% for i=(1:length(delta_ts))
%     delta_t=delta_ts(i);
%     fprintf('delta_t=%d\n',delta_t)
%     T_range=6*3600/delta_t;
%     [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
%     pos_i=[zeros(1,T_range);pos_j];
%     d=srt_dis(pos_i,pos_j);
%     beta=srt_beta(d);
%     gamma=0*ones(J_range,J_range);
%     [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E_means=[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
% ref_E_means=ref_Es(1,:)./J_range;
% save('.\experiments\delta_ts_15_240','delta_ts','srt_E_means','ref_E_means')
% plot(delta_ts,srt_E_means(1,:),'r',delta_ts,srt_E_means(2,:),'r--o',delta_ts,ref_E_means(1,:),'g');

%% interval & J_range 变化
% %N=3;
% J_ranges=(8:8:100);
% intervals=(T_range*delta_t-2*3600)./J_ranges*4;
% ratio=[0,1/4,1/4];
% srt_Es=zeros(2,length(intervals));ref_Es=zeros(1,length(intervals));
% srt_unsatisfied=zeros(1,length(intervals));ref_unsatisfied=zeros(1,length(intervals));
% delta_t=60;
% T_range=6*3600/delta_t;
% for i=(1:length(intervals))
%     interval=intervals(i);
%     J_range=J_ranges(i);
%     fprintf('interval=%2.1fmin, J_range=%d\n',interval/60,J_range)
%     [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
%     pos_i=[zeros(1,T_range);pos_j];
%     d=srt_dis(pos_i,pos_j);
%     beta=srt_beta(d);
%     gamma=0*ones(J_range,J_range);
%     [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E_means=[srt_Es(1,:)./J_ranges;srt_Es(2,:)./J_ranges];
% ref_E_means=ref_Es(1,:)./J_ranges;
% save('.\experiments\J_ranges_8_100','intervals','J_ranges','srt_E_means','ref_E_means','srt_unsatisfied','ref_unsatisfied')
% plot(J_ranges(~srt_unsatisfied),srt_E_means(1,~srt_unsatisfied),'r',J_ranges(~srt_unsatisfied),srt_E_means(2,~srt_unsatisfied),'r--o')%,J_ranges(~ref_unsatisfied),ref_E_means(1,~ref_unsatisfied),'g');

%% N 变化
% Ns=[1,2,3,5,8,10];
% srt_Es=zeros(2,length(Ns));ref_Es=zeros(1,length(Ns));
% srt_unsatisfied=zeros(1,length(Ns));ref_unsatisfied=zeros(1,length(Ns));
% for i=(1:length(Ns))
%     N=Ns(i);
%     fprintf('N=%d\n',N)
%     [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E_means=[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
% ref_E_means=ref_Es./J_range;
% save('.\experiments\Ns_1_10','Ns','srt_E_means','ref_E_means')
% plot(Ns,srt_E_means(1,:),'r',Ns,srt_E_means(2,:),'r--o',Ns,ref_E_means,'g');

%% T_range 变化
T_ranges=[1,(15:15:360)];
srt_Es=zeros(2,length(T_ranges));ref_Es=zeros(1,length(T_ranges));
srt_unsatisfied=zeros(1,length(T_ranges));ref_unsatisfied=zeros(1,length(T_ranges));
C_qos=C_qos/10;
[ref_Es,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,gamma,delta_t,BW);
for i=(1:length(T_ranges))
    T_range_partial=T_ranges(i);
    div=T_range/T_range_partial;
    C_unsatisfied=C_qos;
    fprintf('T_range=%d\n',T_range_partial)
    J_range=(T_range*delta_t-2*3600)/interval*2;
    ii=1;
    while ii<=div && sum(C_unsatisfied)
        pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
        pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
        beta_partial=beta(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
        gamma=0*ones(J_range,J_range);
        srt_Es_partial=zeros(2,1);
        [srt_Es_partial(1),srt_Es_partial(2),C_unsatisfied]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied,N,gamma,delta_t,BW);
        srt_Es(:,i)=srt_Es(:,i)+srt_Es_partial;
        ii=ii+1;
    end
end
srt_E_means=[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
ref_E_means=ref_Es*ones(1,length(T_ranges))./J_range;
save('.\experiments\T_ranges_1_240','T_ranges','srt_E_means','ref_E_means')
plot(T_ranges,srt_E_means(1,:),'r',T_ranges,srt_E_means(2,:),'r--o',T_ranges,ref_E_means,'g');

%% C_qos 变化
% multis=[0.1,0.2,0.5,0.8,1,1.2,1.5,2,2.1,2.2,2.3,2.4,2.5,2.55,2.6,2.7,2.8,2.9,3];
% srt_Es=zeros(2,length(multis));ref_Es=zeros(1,length(multis));
% srt_unsatisfied=zeros(1,length(multis));ref_unsatisfied=zeros(1,length(multis));
% for i=(1:length(multis))
%     multi=multis(i);
%     C_qos_new=C_qos*multi;
%     fprintf('multi=%2.2f\n',multi)
%     [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,gamma,delta_t,BW);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,gamma,delta_t,BW);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E_means=[srt_Es(1,:)./J_range./(multis);srt_Es(2,:)./J_range./(multis)];
% ref_E_means=ref_Es./J_range./(multis);
% save('.\experiments\C_qos_0.1_3.mat','multis','srt_E_means','ref_E_means','srt_unsatisfied','ref_unsatisfied')
% plot(multis(~srt_unsatisfied),srt_E_means(1,~srt_unsatisfied),'r',multis(~srt_unsatisfied),srt_E_means(2,~srt_unsatisfied),'r--o',multis(~ref_unsatisfied),ref_E_means(~ref_unsatisfied),'g');
