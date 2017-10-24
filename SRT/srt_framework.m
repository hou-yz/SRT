%%
clear
clc
interval=15*60;
ratio=[0,1/2,1/2];
sigma2=10^-14;
delta_t=60;
T_range=6*3600/delta_t;%an hour
J_range=(T_range*delta_t-2*3600)/interval*2;
h0_constant=ones(J_range+1,J_range,T_range);
h0_multi=100;
h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
N=3;
P_i_max=10;P_j_max=1;
BW=2*10^6;%bandwidth=2MHz

[pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
pos_i=[zeros(1,T_range);pos_j];
d=srt_dis(pos_i,pos_j);
beta=srt_beta(d);

% service begin time -- bt
% service end time   -- et
bt=ones(1,J_range);et=ones(1,J_range)*T_range;
for j=(1:J_range)
    bt(j)=find(reshape(beta(1,j,:),[1,T_range]),1);
    et(j)=find(reshape(beta(1,j,:),[1,T_range]),1,'last');
end

%%
% 0: h0һά
% 1: h0=ones(J_range+1,J_range,T_range)
% s: h0~CN(0,I)

%%
% N=1;ratio=[0,1/4,1/4];[pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
% pos_i=[zeros(1,T_range);pos_j];
% d=srt_dis(pos_i,pos_j);
% % d1=reshape(d(1,:,:),J_range,T_range);
% beta=srt_beta(d);
% 
% % service begin time -- bt
% % service end time   -- et
% bt=zeros(1,J_range);et=zeros(1,J_range);
% for j=(1:J_range)
%     bt(j)=find(reshape(beta(1,j,:),[1,T_range]),1);
%     et(j)=find(reshape(beta(1,j,:),[1,T_range]),1,'last');
% end
% [srt_E_1,srt_E_2,C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
% srt_unsatisfied=(sum(C_unsatisfied)>0);
% [ref_Es,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
% ref_unsatisfied=(sum(C_unsatisfied)>0);

%% C_qos 
multis=[1,(2:2:20),23,26,30];
srt_E0=zeros(2,length(multis));ref_E0=zeros(1,length(multis));
srt_E1=zeros(2,length(multis));ref_E1=zeros(1,length(multis));
srt_Es=zeros(2,length(multis));ref_Es=zeros(1,length(multis));

srt_unsatisfied=zeros(1,length(multis));ref_unsatisfied=zeros(1,length(multis));
for i=(1:length(multis))
    multi=multis(i);
    C_qos_new=C_qos*multi;
    fprintf('multi=%2.2f\n',multi)
    [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,1,bt,et);
    [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,1,bt,et);
    [srt_E1(1,i),srt_E1(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0_constant,bt,et);
    [ref_E1(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0_constant,bt,et);
    
    [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0,bt,et);
    srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
    [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0,bt,et);
    ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
end
srt_E0_means=[srt_E0(1,:)./J_range./(multis);srt_E0(2,:)./J_range./(multis)];
ref_E0_means=ref_E0./J_range./(multis);
srt_E1_means=[srt_E1(1,:)./J_range./(multis);srt_E1(2,:)./J_range./(multis)];
ref_E1_means=ref_E1./J_range./(multis);
srt_Es_means=[srt_Es(1,:)./J_range./(multis);srt_Es(2,:)./J_range./(multis)];
ref_Es_means=ref_Es./J_range./(multis);
save('.\experiments\C_qos_0.01_50.mat','multis','srt_E0_means','ref_E0_means','srt_E1_means','ref_E1_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
plot(multis(~srt_unsatisfied),srt_E0_means(1,~srt_unsatisfied),'r:^',multis(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied),'b:^',multis(~ref_unsatisfied),ref_E0_means(~ref_unsatisfied),'g:^')
hold on
plot(multis(~srt_unsatisfied),srt_E1_means(1,~srt_unsatisfied),'r--',multis(~srt_unsatisfied),srt_E1_means(2,~srt_unsatisfied),'b--',multis(~ref_unsatisfied),ref_E1_means(~ref_unsatisfied),'g--')
plot(multis(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r-',multis(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'b-',multis(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g-');
hold off
xlabel('C_Q_o_S(Gbit)')
ylabel('Energy Consumption per User per Gbit (J)')
legend('estimation in (2c), cellular-only','estimation in (2c), D2D underlaid','estimation in (2c), reference','h0 constant, cellular-only','h0 constant, D2D underlaid','h0 constant, reference','complete channel, cellular-only','complete channel, D2D underlaid','complete channel, reference')

% %% N 
% Ns=[1,2,3,4,5,6,8,10,12,15,20];
% srt_E0_means=zeros(2,length(Ns));ref_E0_means=zeros(1,length(Ns));
% srt_E1_means=zeros(2,length(Ns));ref_E1_means=zeros(1,length(Ns));
% srt_Es_means=zeros(2,length(Ns));ref_Es_means=zeros(1,length(Ns));
% for ii=(1:3)
%     J_range=(T_range*delta_t-2*3600)./interval*4;
%     h0_constant=ones(J_range+1,J_range,T_range);
%     h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
%     ratio=[0,1/4,1/4];
%     delta_t=4*60;
%     T_range=6*3600/delta_t;%an hour
%     [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
%     C_qos=C_qos*10;
%     pos_i=[zeros(1,T_range);pos_j];
%     d=srt_dis(pos_i,pos_j);
%     beta=srt_beta(d);
%     bt=zeros(1,J_range);et=zeros(1,J_range);
%     for j=(1:J_range)
%         bt(j)=find(reshape(beta(1,j,:),[1,T_range]),1);
%         et(j)=find(reshape(beta(1,j,:),[1,T_range]),1,'last');
%     end
% 
%     srt_E0=zeros(2,length(Ns));ref_E0=zeros(1,length(Ns));
%     srt_E1=zeros(2,length(Ns));ref_E1=zeros(1,length(Ns));
%     srt_Es=zeros(2,length(Ns));ref_Es=zeros(1,length(Ns));
%     srt_unsatisfied=zeros(1,length(Ns));ref_unsatisfied=zeros(1,length(Ns));
%     for i=(1:length(Ns))
%         N=Ns(i);
%         fprintf('N=%d\n',N)
%         [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
%         [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
%         [srt_E1(1,i),srt_E1(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
%         [ref_E1(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
%         
%         [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
%         srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%         [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
%         ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     end
%     srt_E0_means=srt_E0_means+[srt_E0(1,:)./J_range;srt_E0(2,:)./J_range];
%     ref_E0_means=ref_E0_means+ref_E0./J_range;
%     srt_E1_means=srt_E1_means+[srt_E1(1,:)./J_range;srt_E1(2,:)./J_range];
%     ref_E1_means=ref_E1_means+ref_E1./J_range;
%     srt_Es_means=srt_Es_means+[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
%     ref_Es_means=ref_Es_means+ref_Es./J_range;
% end
% srt_E0_means=srt_E0_means./ii;
% ref_E0_means=ref_E0_means./ii;
% srt_E1_means=srt_E1_means./ii;
% ref_E1_means=ref_E1_means./ii;
% srt_Es_means=srt_Es_means./ii;
% ref_Es_means=ref_Es_means./ii;
% save('.\experiments\Ns_1_20','Ns','srt_E0_means','ref_E0_means','srt_E1_means','ref_E1_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
% plot(Ns(~srt_unsatisfied),srt_E0_means(1,~srt_unsatisfied),'r:^',Ns(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied),'b:^',Ns(~ref_unsatisfied),ref_E0_means(~ref_unsatisfied),'g:^')
% hold on
% plot(Ns(~srt_unsatisfied),srt_E1_means(1,~srt_unsatisfied),'r--',Ns(~srt_unsatisfied),srt_E1_means(2,~srt_unsatisfied),'b--',Ns(~ref_unsatisfied),ref_E1_means(~ref_unsatisfied),'g--')
% plot(Ns(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r-',Ns(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'b-',Ns(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g-');
% hold off
% legend('estimation in (2c), cellular-only','estimation in (2c), D2D underlaid','estimation in (2c), reference','h0 constant, cellular-only','h0 constant, D2D underlaid','h0 constant, reference','complete channel, cellular-only','complete channel, D2D underlaid','complete channel, reference')
% xlabel('Number of Subcarriers')
% ylabel('Energy Consumption per User (J)')
% 
% %% T_range
% delta_t=60;
% T_range=6*3600/delta_t;
% [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
% pos_i=[zeros(1,T_range);pos_j];
% d=srt_dis(pos_i,pos_j);
% beta=srt_beta(d);
% bt=zeros(1,J_range);et=zeros(1,J_range);
% for j=(1:J_range)
%     bt(j)=find(reshape(beta(1,j,:),[1,T_range]),1);
%     et(j)=find(reshape(beta(1,j,:),[1,T_range]),1,'last');
% end
% 
% T_ranges=[1,(15:15:T_range)];
% srt_Es=zeros(2,length(T_ranges));srt_E0=zeros(2,length(T_ranges));srt_E1=zeros(2,length(T_ranges));
% srt_unsatisfied=zeros(1,length(T_ranges));ref_unsatisfied=zeros(1,length(T_ranges));
% h0_constant=ones(J_range+1,J_range,T_range);
% h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
% [ref_E0,~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
% [ref_E1,~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
% [ref_Es,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
% for j=(1:3)
%     
%     for i=(1:length(T_ranges))
%         T_range_partial=T_ranges(i);
%         div=T_range/T_range_partial;
%         C_unsatisfied0=C_qos;C_unsatisfied1=C_qos;C_unsatisfieds=C_qos;
%         fprintf('T_range=%d\n',T_range_partial)
%         J_range=(T_range*delta_t-2*3600)/interval*2;
%         h0_constant=ones(J_range+1,J_range,T_range_partial);
%         ii=1;
%         while ii<=ceil(div) && (sum(C_unsatisfied0) || sum(C_unsatisfied1) || sum(C_unsatisfieds))
%             T_range_partial=min(T_range_partial,T_range-T_range_partial*(ii-1));
%             h0_constant=ones(J_range+1,J_range,T_range_partial);
%             h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range_partial*h0_multi])+1j*randn([J_range+1,J_range,T_range_partial*h0_multi]));
%             pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             beta_partial=beta(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             
%             bt_partial=zeros(1,J_range);et_partial=zeros(1,J_range);
%             for jj=(1:J_range)
%                 bt_partial(jj)=max(min(bt(jj),T_range_partial*ii),T_range_partial*(ii-1)+1);
%                 et_partial(jj)=max(min(et(jj),T_range_partial*ii),T_range_partial*(ii-1)+1);
%             end
%             bt_partial=bt_partial-T_range_partial*(ii-1);et_partial=et_partial-T_range_partial*(ii-1);
%             et_partial(bt>T_range_partial*ii)=0;
%             
%             srt_E0_partial=zeros(2,1);
%             [srt_E0_partial(1),srt_E0_partial(2),C_unsatisfied0]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied0,N,delta_t,BW,sigma2,1,bt_partial,et_partial);
%             srt_E0(:,i)=srt_E0(:,i)+srt_E0_partial;
%             
%             srt_E1_partial=zeros(2,1);
%             [srt_E1_partial(1),srt_E1_partial(2),C_unsatisfied1]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied1,N,delta_t,BW,sigma2,h0_constant,bt_partial,et_partial);
%             srt_E1(:,i)=srt_E1(:,i)+srt_E1_partial;
%             
%             srt_Es_partial=zeros(2,1);
%             [srt_Es_partial(1),srt_Es_partial(2),C_unsatisfieds]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfieds,N,delta_t,BW,sigma2,h0,bt_partial,et_partial);
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
% srt_E1_means=[srt_E1(1,:)./J_range/j;srt_E1(2,:)./J_range/j];
% ref_E1_means=ref_E1*ones(1,length(T_ranges))./J_range;
% srt_Es_means=[srt_Es(1,:)./J_range/j;srt_Es(2,:)./J_range/j];
% ref_Es_means=ref_Es*ones(1,length(T_ranges))./J_range;
% save('.\experiments\T_ranges_1_360','T_ranges','srt_E0_means','ref_E0_means','srt_E1_means','ref_E1_means','srt_Es_means','ref_Es_means')
% plot(T_ranges/360,srt_E0_means(1,:),'r:^',T_ranges/360,srt_E0_means(2,:),'b:^',T_ranges/360,ref_E0_means,'g:^',T_ranges/360,srt_E1_means(1,:),'r--',T_ranges/360,srt_E1_means(2,:),'b--',T_ranges/360,ref_E1_means,'g--',T_ranges/360,srt_Es_means(1,:),'r-',T_ranges/360,srt_Es_means(2,:),'b-',T_ranges/360,ref_Es_means,'g-');
% xlabel('Ratio of T_a_c_q to Total Service Duration')
% ylabel('Energy Consumption per User (J)')
% legend('estimation in (2c), cellular-only','estimation in (2c), D2D underlaid','estimation in (2c), reference','h0 constant, cellular-only','h0 constant, D2D underlaid','h0 constant, reference','complete channel, cellular-only','complete channel, D2D underlaid','complete channel, reference')
% 
% %% h0
% T_range=2*3600/delta_t;
% J_range=8;
% h0=ones(J_range+1,J_range,T_range);
% [pos_j,C_qos0]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
% pos_i=[zeros(1,T_range);pos_j];
% d=srt_dis(pos_i,pos_j);
% beta=srt_beta(d);
% bt=zeros(1,J_range);et=zeros(1,J_range);
% for j=(1:J_range)
%     bt(j)=find(reshape(beta(1,j,:),[1,T_range]),1);
%     et(j)=find(reshape(beta(1,j,:),[1,T_range]),1,'last');
% end
% 
% h0_multi=1000;%
% sigma2s=10.^linspace(-16,-9,10);
% srt_Es=zeros(2,length(sigma2s));ref_Es=zeros(1,length(sigma2s));
% srt_E1=zeros(2,length(sigma2s));ref_E1=zeros(1,length(sigma2s));
% srt_E0=zeros(2,length(sigma2s));ref_E0=zeros(1,length(sigma2s));
% srt_unsatisfied=zeros(1,length(sigma2s));ref_unsatisfied=zeros(1,length(sigma2s));
% 
% repeat_max=1;
% for i=(1:length(sigma2s))
%     sigma2=sigma2s(i)
%     C_qos=C_qos0/sqrt(sigma2)*sqrt(10^-14);%abs(C_qos0/log(10^-14/sigma2)*log(10));%C_qos0;%
%     if sigma2==10^-14
%         C_qos=C_qos0;
%     end
%     h0_constant=ones(J_range+1,J_range,T_range);
%     [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
%     [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
%     [srt_E1(1,i),srt_E1(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
%     [ref_E1(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
%     
%     srt_Es_partial=zeros(2,repeat_max);ref_Es_partial=zeros(1,repeat_max);
%     for ii=(1:repeat_max)
%         
%         h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
%         [srt_Es_partial(1,ii),srt_Es_partial(2,ii),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
%         srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%         [ref_Es_partial(1,ii),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);    
%         ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% %         rs(i)=C_unsatisfied;
%     end
%     srt_Es(:,i)=mean(srt_Es_partial,2);ref_Es(:,i)=mean(ref_Es_partial,2);
%     srt_E0(:,i)=srt_E0(:,i)/J_range/C_qos(1);ref_E0(:,i)=ref_E0(:,i)/J_range/C_qos(1);
%     srt_E1(:,i)=srt_E1(:,i)/J_range/C_qos(1);ref_E1(:,i)=ref_E1(:,i)/J_range/C_qos(1);
%     srt_Es(:,i)=srt_Es(:,i)/J_range/C_qos(1);ref_Es(:,i)=ref_Es(:,i)/J_range/C_qos(1); 
% end
% 
% 
% save('.\experiments\sigma2s','sigma2s','srt_Es','ref_Es','srt_E0','ref_E0','srt_E1','ref_E1') 
% 
% 
% plot(sigma2s,srt_E0(1,~srt_unsatisfied),'r:^',sigma2s,srt_E0(2,~srt_unsatisfied),'b:^',sigma2s,ref_E0(~ref_unsatisfied),'g:^')
% hold on
% plot(sigma2s,srt_E1(1,~srt_unsatisfied),'r--',sigma2s,srt_E1(2,~srt_unsatisfied),'b--',sigma2s,ref_E1(~ref_unsatisfied),'g--')
% plot(sigma2s,srt_Es(1,~srt_unsatisfied),'r-',sigma2s,srt_Es(2,~srt_unsatisfied),'b-',sigma2s,ref_Es(~ref_unsatisfied),'g-');
% hold off
% legend('estimation in (2c), cellular-only','estimation in (2c), D2D underlaid','estimation in (2c), reference','h0 constant, cellular-only','h0 constant, D2D underlaid','h0 constant, reference','complete channel, cellular-only','complete channel, D2D underlaid','complete channel, reference')
% 
% xlabel('Sigma^2')
% ylabel('Energy Consumption per User per Gbit (J)')
% 
% %% delay tolerant
% delta_t=60*2;
% T_range=6*3600/delta_t;%an hour
% J_range=(T_range*delta_t-2*3600)/interval*2;
% h0_constant=ones(J_range+1,J_range,T_range);
% h0_multi=100;
% h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
% N=3;
% P_i_max=10;P_j_max=1;
% BW=2*10^6;%bandwidth=2MHz
% 
% [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
% pos_i=[zeros(1,T_range);pos_j];
% d=srt_dis(pos_i,pos_j);
% beta=srt_beta(d);
% 
% bt=ones(1,J_range);et=ones(1,J_range)*T_range;
% for j=(1:J_range)
%     bt(j)=find(reshape(beta(1,j,:),[1,T_range]),1);
%     et(j)=find(reshape(beta(1,j,:),[1,T_range]),1,'last');
% end
% 
% 
% bt0=bt;et0=et;service_time=et0-bt0;
% delay_tolerances=[(1:-0.1:0.1),0.05,0.01];
% 
% srt_E0_means=zeros(2,length(delay_tolerances));ref_E0_means=zeros(1,length(delay_tolerances));
% srt_E1_means=zeros(2,length(delay_tolerances));ref_E1_means=zeros(1,length(delay_tolerances));
% srt_Es_means=zeros(2,length(delay_tolerances));ref_Es_means=zeros(1,length(delay_tolerances));
% 
% ii_max=10;
% for ii=(1:ii_max)
%     h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
%     srt_E0=zeros(2,length(delay_tolerances));ref_E0=zeros(1,length(delay_tolerances));
%     srt_E1=zeros(2,length(delay_tolerances));ref_E1=zeros(1,length(delay_tolerances));
%     srt_Es=zeros(2,length(delay_tolerances));ref_Es=zeros(1,length(delay_tolerances));
%     srt_unsatisfied=zeros(1,length(delay_tolerances));ref_unsatisfied=zeros(1,length(delay_tolerances));
%     for i=(1:length(delay_tolerances))
%         delay_tolerance=delay_tolerances(i);
%         fprintf('delay_tolerance=%2.2f\n',delay_tolerance)
%         
%         bt=round(ii/ii_max*(1-delay_tolerance)*service_time)+bt0;%round(random('uniform',0,(1-delay_tolerance)*service_time))+bt0;%
%         et=bt+round(service_time*delay_tolerance);
%         
%         [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
%         [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,1,bt,et);
%         [srt_E1(1,i),srt_E1(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
%         [ref_E1(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant,bt,et);
%         
%         [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
%         srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%         [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et);
%         ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     end
%     srt_E0_means=srt_E0_means+[srt_E0(1,:)./J_range;srt_E0(2,:)./J_range];
%     ref_E0_means=ref_E0_means+ref_E0./J_range;
%     srt_E1_means=srt_E1_means+[srt_E1(1,:)./J_range;srt_E1(2,:)./J_range];
%     ref_E1_means=ref_E1_means+ref_E1./J_range;
%     srt_Es_means=srt_Es_means+[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
%     ref_Es_means=ref_Es_means+ref_Es./J_range;
% end
% srt_E0_means=srt_E0_means./ii;
% ref_E0_means=ref_E0_means./ii;
% srt_E1_means=srt_E1_means./ii;
% ref_E1_means=ref_E1_means./ii;
% srt_Es_means=srt_Es_means./ii;
% ref_Es_means=ref_Es_means./ii;
% save('.\experiments\delays','delay_tolerance','srt_E0_means','ref_E0_means','srt_E1_means','ref_E1_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
% plot(delay_tolerances(~srt_unsatisfied),srt_E0_means(1,~srt_unsatisfied),'r:^',delay_tolerances(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied),'b:^')%,delay_tolerances(~ref_unsatisfied),ref_E0_means(~ref_unsatisfied),'g:^')
% hold on
% plot(delay_tolerances(~srt_unsatisfied),srt_E1_means(1,~srt_unsatisfied),'r--',delay_tolerances(~srt_unsatisfied),srt_E1_means(2,~srt_unsatisfied),'b--')%,delay_tolerances(~ref_unsatisfied),ref_E1_means(~ref_unsatisfied),'g--')
% plot(delay_tolerances(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r-',delay_tolerances(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'b-')%,delay_tolerances(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g-');
% hold off
% legend('estimation in (2c), cellular-only','estimation in (2c), D2D underlaid','h0 constant, cellular-only','h0 constant, D2D underlaid','complete channel, cellular-only','complete channel, D2D underlaid')
% xlabel('Ratio of Delay Tolerance to Service Time')
% ylabel('Energy Consumption per User (J)')
% 
