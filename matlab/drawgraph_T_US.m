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

%% T_range



%% 编辑这里改横轴分辨率
T_ranges=linspace(1,T_range);
%%



srt_Es=zeros(2,length(T_ranges));srt_E0=zeros(2,length(T_ranges));
srt_unsatisfied=zeros(1,length(T_ranges));ref_unsatisfied=zeros(1,length(T_ranges));
[ref_E0,~]=ref_algorithm(beta0,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt0,et0);
[ref_Es,C_unsatisfied]=ref_algorithm(betaS,J_range,T_range*t_step,P_i_max,P_j_max,C_qos,N,delta_t/t_step,BW,sigma2,h0,btS,etS);
for j=(1:10)
    
    for i=(1:length(T_ranges))
        T_range_partial=T_ranges(i);
        div=T_range/T_range_partial;
        C_unsatisfied0=C_qos;C_unsatisfieds=C_qos;
        fprintf('T_range=%d\n',T_range_partial)
        
        ii=1;
        while ii<=ceil(div) && (sum(C_unsatisfied0) || sum(C_unsatisfieds))
            T_range_partial=min(T_range_partial,T_range-T_range_partial*(ii-1));
            pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
            pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
            beta0_partial=beta0(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
            betaS_partial=betaS(:,:,T_range_partial*t_step*(ii-1)+1 : T_range_partial*t_step*ii);
            h0_partial=h0(:,:,T_range_partial*t_step*(ii-1)*10+1 : T_range_partial*t_step*10*ii);
            
            bt_partial=zeros(1,J_range);et_partial=zeros(1,J_range);
            for jj=(1:J_range)
                bt_partial(jj)=max(min(bt0(jj),T_range_partial*ii),T_range_partial*(ii-1)+1);
                et_partial(jj)=max(min(et0(jj),T_range_partial*ii),T_range_partial*(ii-1)+1);
            end
            bt_partial=bt_partial-T_range_partial*(ii-1);et_partial=et_partial-T_range_partial*(ii-1);
            et_partial(bt0>T_range_partial*ii)=0;
            
            srt_E0_partial=zeros(2,1);
            [~,srt_E0_partial(2),C_unsatisfied0]=srt_algorithm_v2(beta0_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied0,N,delta_t,BW,sigma2,1,bt_partial,et_partial);
            srt_E0(:,i)=srt_E0(:,i)+srt_E0_partial;
            
            bt_partial=zeros(1,J_range);et_partial=zeros(1,J_range);
            for jj=(1:J_range)
                bt_partial(jj)=max(min(btS(jj),T_range_partial*t_step*ii),T_range_partial*(ii-1)*t_step+1);
                et_partial(jj)=max(min(etS(jj),T_range_partial*t_step*ii),T_range_partial*(ii-1)*t_step+1);
            end
            bt_partial=bt_partial-T_range_partial*(ii-1)*t_step;et_partial=et_partial-T_range_partial*(ii-1)*t_step;
            et_partial(bt0>T_range_partial*ii)=0;
            srt_Es_partial=zeros(2,1);
            [~,srt_Es_partial(2),C_unsatisfieds]=srt_algorithm_v2(betaS_partial,J_range,T_range_partial*t_step,P_i_max,P_j_max,C_unsatisfieds,N,delta_t/t_step,BW,sigma2,h0_partial,bt_partial,et_partial);
            srt_Es(:,i)=srt_Es(:,i)+srt_Es_partial;
            
            ii=ii+1;
        end
        
    end
    
end

srt_E0_means=[srt_E0(1,:)./J_range/j;srt_E0(2,:)./J_range/j];
ref_E0_means=ref_E0*ones(1,length(T_ranges))./J_range;
srt_Es_means=[srt_Es(1,:)./J_range/j;srt_Es(2,:)./J_range/j];
ref_Es_means=ref_Es*ones(1,length(T_ranges))./J_range;
save('.\experiments\T_ranges_1_360','T_ranges','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means')
i=[1,3,5,7,9,11];
plot(T_ranges(i)/120,srt_E0_means(2,i)/1000,'r-^',T_ranges(i)/120,srt_Es_means(2,i)/1000,'b-o',T_ranges(i)/120,ref_Es_means(i)/1000,'g-*');
grid on
xlabel('Ratio of T^U^S to Total Service Duration')
ylabel('Energy Consumption per User (kJ)')
legend('large-scale CSI, ship-to-ship/shore','(genius-aided) full CSI, ship-to-ship/shore','without CSI, round-robin')
