%%
clear
clc
interval=15*60;
ratio=[0,1/2,1/2];
sigma2=10^-10.7;
delta_t=60;
T_range=2*3600/delta_t;%6*3600/delta_t;
J_range=8;%(T_range*delta_t-2*3600)/interval*2;
t_step=1;
h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*t_step*10])+1j*randn([J_range+1,J_range,T_range*t_step*10]));
N=32;
P_i_max=10;P_j_max=1;
BW=2*10^6;%bandwidth=2MHz

[pos_j,~]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
pos_i=[zeros(1,T_range);pos_j];
d=srt_dis(pos_i,pos_j);
beta0=srt_beta(d);
bt0=ones(1,J_range);et0=ones(1,J_range)*T_range;
for j=(1:J_range)
    bt0(j)=find(reshape(beta0(1,j,:),[1,T_range]),1);
    et0(j)=find(reshape(beta0(1,j,:),[1,T_range]),1,'last');
end
[pos_j,C_qos]=srt_pos_and_qos(J_range,T_range*t_step,delta_t/t_step,interval,ratio);
pos_i=[zeros(1,T_range*t_step);pos_j];
d=srt_dis(pos_i,pos_j);
betaS=srt_beta(d);
btS=ones(1,J_range);etS=ones(1,J_range)*T_range;
for j=(1:J_range)
    btS(j)=find(reshape(betaS(1,j,:),[1,T_range*t_step]),1);
    etS(j)=find(reshape(betaS(1,j,:),[1,T_range*t_step]),1,'last');
end

% service begin time -- bt0
% service end time   -- et0

%%
% 0: h0һά
% 1: h0=ones(J_range+1,J_range,T_range)
% s: h0~CN(0,I)


%% C_qos 
% multis=linspace(0.2,3,20);
% srt_E0=zeros(2,length(multis));ref_E0=zeros(1,length(multis));
% srt_Es=zeros(2,length(multis));ref_Es=zeros(1,length(multis));
% 
% srt_unsatisfied=zeros(1,length(multis));ref_unsatisfied=zeros(1,length(multis));
% for i=(1:length(multis))
%     multi=multis(i);
%     C_qos_new=C_qos*multi;
%     fprintf('multi=%2.2f\n',multi)
%     [S0,srt_E0(2,i),~]=srt_algorithm_v2(beta0,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,1,bt0,et0);
%     [ref_E0(1,i),~]=ref_algorithm(beta0,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,1,bt0,et0);
% 
%     [Ss,srt_Es(2,i),C_unsatisfied]=srt_algorithm_v2(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos_new,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_Es(1,i),C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos_new,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E0_means=[srt_E0(1,:)./J_range./(multis);srt_E0(2,:)./J_range./(multis)];
% ref_E0_means=ref_E0./J_range./(multis);
% srt_Es_means=[srt_Es(1,:)./J_range./(multis);srt_Es(2,:)./J_range./(multis)];
% ref_Es_means=ref_Es./J_range./(multis);
% % save('.\experiments\C_qos_0.01_50.mat','multis','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
% plot(multis(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied)/1000,'r-*')
% hold on
% plot(multis(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied)/1000,'b-',multis(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied)/1000,'g-');
% hold off
% xlabel('C_Q_o_S(Gbit)')
% ylabel('Energy Consumption per User per Gbit (kJ)')
% legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore','without CSI, round-robin')

%% T_range
% T_ranges=[1,(12:12:T_range)];
% srt_Es=zeros(2,length(T_ranges));srt_E0=zeros(2,length(T_ranges));
% srt_unsatisfied=zeros(1,length(T_ranges));ref_unsatisfied=zeros(1,length(T_ranges));
% [ref_E0,~]=ref_algorithm(beta0,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt0,et0);
% [ref_Es,C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
% for j=(1:4)
%     
%     for i=(1:length(T_ranges))
%         T_range_partial=T_ranges(i);
%         div=T_range/T_range_partial;
%         C_unsatisfied0=C_qos;C_unsatisfieds=C_qos;
%         fprintf('T_range=%d\n',T_range_partial)
%         
%         ii=1;
%         while ii<=ceil(div) && (sum(C_unsatisfied0) || sum(C_unsatisfieds))
%             T_range_partial=min(T_range_partial,T_range-T_range_partial*(ii-1));
%             pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             beta0_partial=beta0(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             betaS_partial=betaS(:,:,T_range_partial*t_step*(ii-1)+1 : T_range_partial*t_step*ii);
%             h0_partial=h0(:,:,T_range_partial*t_step*(ii-1)*10+1 : T_range_partial*t_step*10*ii);
%             
%             bt_partial=zeros(1,J_range);et_partial=zeros(1,J_range);
%             for jj=(1:J_range)
%                 bt_partial(jj)=max(min(bt0(jj),T_range_partial*ii),T_range_partial*(ii-1)+1);
%                 et_partial(jj)=max(min(et0(jj),T_range_partial*ii),T_range_partial*(ii-1)+1);
%             end
%             bt_partial=bt_partial-T_range_partial*(ii-1);et_partial=et_partial-T_range_partial*(ii-1);
%             et_partial(bt0>T_range_partial*ii)=0;
%             
%             srt_E0_partial=zeros(2,1);
%             [~,srt_E0_partial(2),C_unsatisfied0]=srt_algorithm_v2(beta0_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied0,N,delta_t,BW,sigma2,1,bt_partial,et_partial);
%             srt_E0(:,i)=srt_E0(:,i)+srt_E0_partial;
%             
%             bt_partial=zeros(1,J_range);et_partial=zeros(1,J_range);
%             for jj=(1:J_range)
%                 bt_partial(jj)=max(min(btS(jj),T_range_partial*t_step*ii),T_range_partial*(ii-1)*t_step+1);
%                 et_partial(jj)=max(min(etS(jj),T_range_partial*t_step*ii),T_range_partial*(ii-1)*t_step+1);
%             end
%             bt_partial=bt_partial-T_range_partial*(ii-1)*t_step;et_partial=et_partial-T_range_partial*(ii-1)*t_step;
%             et_partial(bt0>T_range_partial*ii)=0;
%             srt_Es_partial=zeros(2,1);
%             [~,srt_Es_partial(2),C_unsatisfieds]=srt_algorithm_v2(betaS_partial,J_range,T_range_partial*t_step,P_i_max,P_j_max,C_unsatisfieds,N,delta_t/t_step,BW,sigma2,h0_partial,bt_partial,et_partial);
%             srt_Es(:,i)=srt_Es(:,i)+srt_Es_partial;
%             
%             ii=ii+1;
%         end
%         
%     end
%     
% end
% 
% srt_E0_means=[srt_E0(1,:)./J_range/j;srt_E0(2,:)./J_range/j];
% ref_E0_means=ref_E0*ones(1,length(T_ranges))./J_range;
% srt_Es_means=[srt_Es(1,:)./J_range/j;srt_Es(2,:)./J_range/j];
% ref_Es_means=ref_Es*ones(1,length(T_ranges))./J_range;
% % save('.\experiments\T_ranges_1_360','T_ranges','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means')
% plot(T_ranges/120,srt_E0_means(2,:)/1000,'r-*',T_ranges/120,srt_Es_means(2,:)/1000,'b-',T_ranges/120,ref_Es_means/1000,'g-');
% xlabel('Ratio of T_a_c_q to Total Service Duration')
% ylabel('Energy Consumption per User (kJ)')
% legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore','without CSI, round-robin')

%% SNR
% C_qos0=C_qos;
% P_i_s=10.^(0.3:0.1:1.5);
% srt_Es=zeros(2,length(P_i_s));ref_Es=zeros(1,length(P_i_s));
% srt_E0=zeros(2,length(P_i_s));ref_E0=zeros(1,length(P_i_s));
% srt_unsatisfied=zeros(1,length(P_i_s));ref_unsatisfied=zeros(1,length(P_i_s));
% 
% repeat_max=1;
% for i=(1:length(P_i_s))
%     P_i=P_i_s(i)
%     P_j=P_i/10;
%     C_qos=C_qos0;%*sqrt(P_i)/sqrt(P_i_max);%abs(C_qos0/log(10^-14/sigma2)*log(10));%
%     [~,srt_E0(2,i),~]=srt_algorithm_v2(beta0,J_range,T_range,P_i,P_j,C_qos,N,delta_t,BW,sigma2,1,bt0,et0);
%     [ref_E0(1,i),~]=ref_algorithm(beta0,J_range,T_range,P_i,P_j,C_qos,N,delta_t,BW,sigma2,1,bt0,et0);
%     srt_Es_partial=zeros(2,repeat_max);ref_Es_partial=zeros(1,repeat_max);
%     for ii=(1:repeat_max)
%         [~,srt_Es_partial(2,ii),C_unsatisfied]=srt_algorithm_v2(betaS,J_range,T_range*t_step,P_i,P_j,C_qos,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
%         srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%         [ref_Es_partial(1,ii),C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i,P_j,C_qos,N,delta_t/t_step,BW,sigma2,h0,btS,etS);    
%         ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     end
%     srt_Es(:,i)=mean(srt_Es_partial,2);ref_Es(:,i)=mean(ref_Es_partial,2);
%     srt_E0(:,i)=srt_E0(:,i)/J_range/C_qos(1)*10^9;ref_E0(:,i)=ref_E0(:,i)/J_range/C_qos(1)*10^9;
%     srt_Es(:,i)=srt_Es(:,i)/J_range/C_qos(1)*10^9;ref_Es(:,i)=ref_Es(:,i)/J_range/C_qos(1)*10^9; 
% end
% 
% % save('.\experiments\SNRs','P_i_s','srt_Es','ref_Es','srt_E0','ref_E0') 
plot(P_i_s,srt_E0(2,:)/1000,'r-*')
hold on
plot(P_i_s,srt_Es(2,:)/1000,'b-')%,P_i_s,ref_Es(:)/1000,'g-');
hold off
legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore')%,'without CSI, round-robin')
xlabel('Transmission Power of BS on given subcarrier (W)')
ylabel('Energy Consumption per User (kJ)')

%% delay tolerant

service_time=et0-bt0;
delay_tolerances=[(1:-0.1:0.1),0.05,0.01];

srt_E0_means=zeros(2,length(delay_tolerances));ref_E0_means=zeros(1,length(delay_tolerances));
srt_Es_means=zeros(2,length(delay_tolerances));ref_Es_means=zeros(1,length(delay_tolerances));

ii_max=1;
for ii=(1:ii_max)
    srt_E0=zeros(2,length(delay_tolerances));ref_E0=zeros(1,length(delay_tolerances));
    srt_Es=zeros(2,length(delay_tolerances));ref_Es=zeros(1,length(delay_tolerances));
    srt_unsatisfied=zeros(1,length(delay_tolerances));ref_unsatisfied=zeros(1,length(delay_tolerances));
    for i=(1:length(delay_tolerances))
        delay_tolerance=delay_tolerances(i);
        fprintf('delay_tolerance=%2.2f\n',delay_tolerance)
        
        bt=round(ii/ii_max*(1-delay_tolerance)*service_time)+bt0;%round(random('uniform',0,(1-delay_tolerance)*service_time))+bt0;%
        et=bt+round(service_time*delay_tolerance);
        [~,srt_E0(2,i),~]=srt_algorithm_v2(beta0,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
        [ref_E0(1,i),~]=ref_algorithm(beta0,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
        
        bt=round(ii/ii_max*(1-delay_tolerance)*service_time)*t_step+btS;%round(random('uniform',0,(1-delay_tolerance)*service_time))+bt0;%
        et=bt+round(service_time*delay_tolerance)*t_step;
        [~,srt_Es(2,i),C_unsatisfied]=srt_algorithm_v2(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos,N,delta_t/t_step,BW,sigma2,h0,bt,et);
        srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
        [ref_Es(1,i),C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos,N,delta_t/t_step,BW,sigma2,h0,bt,et);
        ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
    end
    srt_E0_means=srt_E0_means+[srt_E0(1,:)./J_range;srt_E0(2,:)./J_range];
    ref_E0_means=ref_E0_means+ref_E0./J_range;
    srt_Es_means=srt_Es_means+[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
    ref_Es_means=ref_Es_means+ref_Es./J_range;
end
srt_E0_means=srt_E0_means./ii;
ref_E0_means=ref_E0_means./ii;
srt_Es_means=srt_Es_means./ii;
ref_Es_means=ref_Es_means./ii;
% save('.\experiments\delays','delay_tolerances','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
plot(delay_tolerances(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied)/1000,'r-*')%,delay_tolerances(~ref_unsatisfied),ref_E0_means(~ref_unsatisfied),'g:*')
hold on
plot(delay_tolerances(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied)/1000,'b-')%,delay_tolerances(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g-');
hold off
legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore')
xlabel('Ratio of Delay Tolerance to Service Time')
ylabel('Energy Consumption per User (kJ)')

