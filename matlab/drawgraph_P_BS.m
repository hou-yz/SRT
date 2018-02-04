%%
clc
clear
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
% s: h0~CN(0,I)


%% SNR
C_qos0=C_qos;
P_i_s=10.^(0.3:0.1:1.5);
srt_Es=zeros(2,length(P_i_s));ref_Es=zeros(1,length(P_i_s));
srt_E0=zeros(2,length(P_i_s));ref_E0=zeros(1,length(P_i_s));
srt_unsatisfied=zeros(1,length(P_i_s));ref_unsatisfied=zeros(1,length(P_i_s));

repeat_max=1;
for i=(1:length(P_i_s))
    P_i=P_i_s(i)
    P_j=P_i/10;
    C_qos=C_qos0;%*sqrt(P_i)/sqrt(P_i_max);%abs(C_qos0/log(10^-14/sigma2)*log(10));%
    [~,srt_E0(2,i),~]=srt_algorithm_v2(beta0,J_range,T_range,P_i,P_j,C_qos,N,delta_t,BW,sigma2,1,bt0,et0);
    [ref_E0(1,i),~]=ref_algorithm(beta0,J_range,T_range,P_i,P_j,C_qos,N,delta_t,BW,sigma2,1,bt0,et0);
    srt_Es_partial=zeros(2,repeat_max);ref_Es_partial=zeros(1,repeat_max);
    for ii=(1:repeat_max)
        [~,srt_Es_partial(2,ii),C_unsatisfied]=srt_algorithm_v2(betaS,J_range,T_range*t_step,P_i,P_j,C_qos,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
        srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
        [ref_Es_partial(1,ii),C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i,P_j,C_qos,N,delta_t/t_step,BW,sigma2,h0,btS,etS);    
        ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
    end
    srt_Es(:,i)=mean(srt_Es_partial,2);ref_Es(:,i)=mean(ref_Es_partial,2);
    srt_E0(:,i)=srt_E0(:,i)/J_range/C_qos(1)*10^9;ref_E0(:,i)=ref_E0(:,i)/J_range/C_qos(1)*10^9;
    srt_Es(:,i)=srt_Es(:,i)/J_range/C_qos(1)*10^9;ref_Es(:,i)=ref_Es(:,i)/J_range/C_qos(1)*10^9; 
end

save('.\experiments\SNRs','P_i_s','srt_Es','ref_Es','srt_E0','ref_E0') 
plot(P_i_s,srt_E0(2,:)/1000,'r-^')
grid on
hold on
plot(P_i_s,srt_Es(2,:)/1000,'b-o')%,P_i_s,ref_Es(:)/1000,'g-');
hold off
legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore')%,'without CSI, round-robin')
xlabel('Transmit Power of BS on given subcarrier (W)')
ylabel('Energy Consumption per User (kJ)')
