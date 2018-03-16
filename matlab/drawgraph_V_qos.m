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
% 0: h0Ò»Î¬
% s: h0~CN(0,I)


%% C_qos 
multis=(0.5:0.5:3);
srt_E0=zeros(2,length(multis));ref_E0=zeros(1,length(multis));
srt_Es=zeros(2,length(multis));ref_Es=zeros(1,length(multis));

srt_unsatisfied=zeros(1,length(multis));ref_unsatisfied=zeros(1,length(multis));
for i=(1:length(multis))
    multi=multis(i);
    C_qos_new=C_qos*multi;
    fprintf('multi=%2.2f\n',multi)
    [S0,srt_E0(2,i),~]=srt_algorithm_v2(beta0,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,1,bt0,et0);
    [ref_E0(1,i),~]=ref_algorithm(beta0,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,1,bt0,et0);

    [Ss,srt_Es(2,i),C_unsatisfied]=srt_algorithm_v2(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos_new,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
    srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
    [ref_Es(1,i),C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos_new,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
    ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
end
srt_E0_means=[srt_E0(1,:)./J_range./(multis);srt_E0(2,:)./J_range./(multis)];
ref_E0_means=ref_E0./J_range./(multis);
srt_Es_means=[srt_Es(1,:)./J_range./(multis);srt_Es(2,:)./J_range./(multis)];
ref_Es_means=ref_Es./J_range./(multis);

%±£´æ
%save('.\experiments\C_qos_0.01_50.mat','multis','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
plot(multis(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied)/1000,'r-^')
grid on
hold on
plot(multis(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied)/1000,'b-o',multis(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied)/1000,'g-*');
hold off
xlabel('V^Q^o^S(Gbit)')
ylabel('Energy Consumption per User per Gbit (kJ)')
legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore','without CSI, round-robin')
