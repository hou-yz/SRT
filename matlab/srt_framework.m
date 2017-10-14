%%
clear
clc
%航线1，航线2，航线3 发船间隔
interval=15*60;%单位：s
%航线1，航线2，航线3 占比
ratio=[0,1/2,1/2];
sigma2=10^-14;
delta_t=60;%单位：秒
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
% d1=reshape(d(1,:,:),J_range,T_range);%到基站距离
beta=srt_beta(d);
% beta1=reshape(beta(1,:,:),J_range,T_range);
% 
% [srt_E_1,srt_E_2,C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
% srt_unsatisfied=(sum(C_unsatisfied)>0);
% [ref_Es,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
% ref_unsatisfied=(sum(C_unsatisfied)>0);

%% C_qos 变化
% multis=[1,2,4,10,20,30,40,50,60];
% srt_E0=zeros(2,length(multis));ref_E0=zeros(1,length(multis));
% srt_Es=zeros(2,length(multis));ref_Es=zeros(1,length(multis));
% 
% srt_unsatisfied=zeros(1,length(multis));ref_unsatisfied=zeros(1,length(multis));
% for i=(1:length(multis))
%     multi=multis(i);
%     C_qos_new=C_qos*multi;
%     fprintf('multi=%2.2f\n',multi)
%     [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0_constant);
%     [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0_constant);
%     
%     [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0);
%     srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos_new,N,delta_t,BW,sigma2,h0);
%     ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
% end
% srt_E0_means=[srt_E0(1,:)./J_range./(multis);srt_E0(2,:)./J_range./(multis)];
% ref_E0_means=ref_E0./J_range./(multis);
% srt_Es_means=[srt_Es(1,:)./J_range./(multis);srt_Es(2,:)./J_range./(multis)];
% ref_Es_means=ref_Es./J_range./(multis);
% save('.\experiments\C_qos_0.01_50.mat','multis','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
% plot(multis(~srt_unsatisfied),srt_E0_means(1,~srt_unsatisfied),'r--',multis(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied),'r--o',multis(~ref_unsatisfied),ref_E0_means(~ref_unsatisfied),'g--',multis(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r-',multis(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'r-o',multis(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g-');
% xlabel('C_Q_o_S(Gbit)')
% ylabel('Energy Consumption per User per Gbit (J)')
% legend('h_0 constant, cellular-only','h_0 constant, D2D underlaid','h_0 constant, reference','h_0 gaussian, cellular-only','h_0 gaussian, D2D underlaid','h_0 gaussian, reference')

%% J_range 变化
% J_ranges=(8:8:100);
% N=2;
% intervals=(T_range*delta_t-2*3600)./J_ranges*4;
% ratio=[0,1/4,1/4];
% srt_E0=zeros(2,length(intervals));ref_E0=zeros(1,length(intervals));
% srt_Es=zeros(2,length(intervals));ref_Es=zeros(1,length(intervals));
% srt_unsatisfied=zeros(1,length(intervals));ref_unsatisfied=zeros(1,length(intervals));
% delta_t=5*60;
% T_range=6*3600/delta_t;
% for i=(1:length(intervals))
%     interval=intervals(i);
%     J_range=J_ranges(i);
%     h0_constant=ones(J_range+1,J_range,T_range);
%     h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
%     fprintf('interval=%2.1fmin, J_range=%d\n',interval/60,J_range)
%     srt_E0_tmps=zeros(2,5);ref_E0_tmps=zeros(1,5);
%     srt_Es_tmps=zeros(2,5);ref_Es_tmps=zeros(1,5);
%     for j=(1:5)
%         [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
%         pos_i=[zeros(1,T_range);pos_j];
%         C_qos=C_qos*5;
%         d=srt_dis(pos_i,pos_j);
%         beta=srt_beta(d);
%         [srt_Es_tmps(1,j),srt_Es_tmps(2,j),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant);
%         [ref_Es_tmps(1,j),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant);
%         [srt_Es_tmps(1,j),srt_Es_tmps(2,j),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
%         srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%         [ref_Es_tmps(1,j),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
%         ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     end
%     srt_E0(:,i)=mean(srt_E0_tmps,2);
%     ref_E0(:,i)=mean(ref_E0_tmps,2);
%     srt_Es(:,i)=mean(srt_Es_tmps,2);
%     ref_Es(:,i)=mean(ref_Es_tmps,2);
% end
% srt_E0_means=[srt_E0(1,:)./J_ranges;srt_E0(2,:)./J_ranges];
% ref_E0_means=ref_E0(1,:)./J_ranges;
% srt_Es_means=[srt_Es(1,:)./J_ranges;srt_Es(2,:)./J_ranges];
% ref_Es_means=ref_Es(1,:)./J_ranges;
% save('.\experiments\J_ranges_8_100','intervals','J_ranges','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
% plot(J_ranges(~srt_unsatisfied),srt_E0_means(1,~srt_unsatisfied),'r',J_ranges(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied),'r--o',J_ranges(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r',J_ranges(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'r--o')%,J_ranges(~ref_unsatisfied),ref_E_means(1,~ref_unsatisfied),'g');
% xlabel('Number of Users')
% ylabel('Energy Consumption per User (J)')
% legend('h_0 constant, cellular-only','h_0 constant, D2D underlaid','h_0 gaussian, cellular-only','h_0 gaussian, D2D underlaid')

%% N 变化
% Ns=[1,2,4,5,6,8,10,12,15,20];
% srt_E0_means=zeros(2,length(Ns));ref_E0_means=zeros(1,length(Ns));
% srt_Es_means=zeros(2,length(Ns));ref_Es_means=zeros(1,length(Ns));
% for ii=(1:4)
%     J_range=(T_range*delta_t-2*3600)./interval*4;
%     h0_constant=ones(J_range+1,J_range,T_range);
%     h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
%     ratio=[0,1/4,1/4];
%     delta_t=4*60;%单位：秒
%     T_range=6*3600/delta_t;%an hour
%     [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
%     C_qos=C_qos*10;
%     pos_i=[zeros(1,T_range);pos_j];
%     d=srt_dis(pos_i,pos_j);
%     beta=srt_beta(d);
%     srt_E0=zeros(2,length(Ns));ref_E0=zeros(1,length(Ns));
%     srt_Es=zeros(2,length(Ns));ref_Es=zeros(1,length(Ns));
%     srt_unsatisfied=zeros(1,length(Ns));ref_unsatisfied=zeros(1,length(Ns));
%     for i=(1:length(Ns))
%         N=Ns(i);
%         fprintf('N=%d\n',N)
%         [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant);
%         [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant);
%         
%         [srt_Es(1,i),srt_Es(2,i),C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
%         srt_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%         [ref_Es(1,i),C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
%         ref_unsatisfied(1,i)=(sum(C_unsatisfied)>0);
%     end
%     srt_E0_means=srt_E0_means+[srt_E0(1,:)./J_range;srt_E0(2,:)./J_range];
%     ref_E0_means=ref_E0_means+ref_E0./J_range;
%     srt_Es_means=srt_Es_means+[srt_Es(1,:)./J_range;srt_Es(2,:)./J_range];
%     ref_Es_means=ref_Es_means+ref_Es./J_range;
% end
% srt_E0_means=srt_E0_means./ii;
% ref_E0_means=ref_E0_means./ii;
% srt_Es_means=srt_Es_means./ii;
% ref_Es_means=ref_Es_means./ii;
% save('.\experiments\Ns_1_20','Ns','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means','srt_unsatisfied','ref_unsatisfied')
% plot(Ns(~srt_unsatisfied),srt_E0_means(1,~srt_unsatisfied),'r--',Ns(~srt_unsatisfied),srt_E0_means(2,~srt_unsatisfied),'r--o',Ns(~ref_unsatisfied),ref_E0_means(~ref_unsatisfied),'g--',Ns(~srt_unsatisfied),srt_Es_means(1,~srt_unsatisfied),'r-',Ns(~srt_unsatisfied),srt_Es_means(2,~srt_unsatisfied),'r-o',Ns(~ref_unsatisfied),ref_Es_means(~ref_unsatisfied),'g-');
% xlabel('Number of Subcarriers')
% ylabel('Energy Consumption per User (J)')
% legend('h_0 constant, cellular-only','h_0 constant, D2D underlaid','h_0 constant, reference','h_0 gaussian, cellular-only','h_0 gaussian, D2D underlaid','h_0 gaussian, reference')

%% T_range 变化
% delta_t=60;
% T_range=6*3600/delta_t;
% [pos_j,C_qos]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
% pos_i=[zeros(1,T_range);pos_j];
% d=srt_dis(pos_i,pos_j);
% beta=srt_beta(d);
% 
% T_ranges=[1,(15:15:T_range)];
% srt_Es=zeros(2,length(T_ranges));srt_E0=zeros(2,length(T_ranges));
% srt_unsatisfied=zeros(1,length(T_ranges));ref_unsatisfied=zeros(1,length(T_ranges));
% h0_constant=ones(J_range+1,J_range,T_range);
% h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range*h0_multi])+1j*randn([J_range+1,J_range,T_range*h0_multi]));
% [ref_E0,~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0_constant);
% [ref_Es,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
% for j=(1:4)
%     
%     for i=(1:length(T_ranges))
%         T_range_partial=T_ranges(i);
%         div=T_range/T_range_partial;
%         C_unsatisfied=C_qos;
%         fprintf('T_range=%d\n',T_range_partial)
%         J_range=(T_range*delta_t-2*3600)/interval*2;
%         h0_constant=ones(J_range+1,J_range,T_range_partial);
%         ii=1;
%         while ii<=ceil(div) && sum(C_unsatisfied)
%             T_range_partial=min(T_range_partial,T_range-T_range_partial*(ii-1));
%             h0_constant=ones(J_range+1,J_range,T_range_partial);
%             pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             beta_partial=beta(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             srt_E0_partial=zeros(2,1);
%             [srt_E0_partial(1),srt_E0_partial(2),C_unsatisfied]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied,N,delta_t,BW,sigma2,h0_constant);
%             srt_E0(:,i)=srt_E0(:,i)+srt_E0_partial;
%             ii=ii+1;
%         end
%         T_range_partial=T_ranges(i);
%         C_unsatisfied=C_qos;
%         ii=1;
%         while ii<=ceil(div) && sum(C_unsatisfied)
%             T_range_partial=min(T_range_partial,T_range-T_range_partial*(ii-1));
%             h0=sqrt(1/2)*(randn([J_range+1,J_range,T_range_partial*h0_multi])+1j*randn([J_range+1,J_range,T_range_partial*h0_multi]));
%             pos_j_partial=pos_j(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             pos_i_partial=pos_i(:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             beta_partial=beta(:,:,T_range_partial*(ii-1)+1 : T_range_partial*ii);
%             srt_Es_partial=zeros(2,1);
%             [srt_Es_partial(1),srt_Es_partial(2),C_unsatisfied]=srt_algorithm(beta_partial,J_range,T_range_partial,P_i_max,P_j_max,C_unsatisfied,N,delta_t,BW,sigma2,h0);
%             srt_Es(:,i)=srt_Es(:,i)+srt_Es_partial;
%             ii=ii+1;
%         end
%     end
%     
% end
% 
% srt_E0_means=[srt_E0(1,:)./J_range/j;srt_E0(2,:)./J_range/j];
% ref_E0_means=ref_E0*ones(1,length(T_ranges))./J_range;
% srt_Es_means=[srt_Es(1,:)./J_range/j;srt_Es(2,:)./J_range/j];
% ref_Es_means=ref_Es*ones(1,length(T_ranges))./J_range;
% save('.\experiments\T_ranges_1_360','T_ranges','srt_E0_means','ref_E0_means','srt_Es_means','ref_Es_means')
% plot(T_ranges/360*100,srt_E0_means(1,:),'r--',T_ranges/360*100,srt_E0_means(2,:),'r--o',T_ranges/360*100,ref_E0_means,'g--',T_ranges/360*100,srt_Es_means(1,:),'r',T_ranges/360*100,srt_Es_means(2,:),'r-o',T_ranges/360*100,ref_Es_means,'g');
% xlabel('Percentage of Time Slots with Acquired CSI')
% ylabel('Energy Consumption per User (J)')
% legend('h_0 constant, cellular-only','h_0 constant, D2D underlaid','h_0 constant, reference','h_0 gaussian, cellular-only','h_0 gaussian, D2D underlaid','h_0 gaussian, reference')

%% h0 快衰落影响
T_range=2*3600/delta_t;
J_range=8;
h0=ones(J_range+1,J_range,T_range);
[pos_j,C_qos0]=srt_pos_and_qos(J_range,T_range,delta_t,interval,ratio);
pos_i=[zeros(1,T_range);pos_j];
d=srt_dis(pos_i,pos_j);
beta=srt_beta(d);

h0_multi=1000;%
sigma2s=10.^(-18:-12);
srt_Es=zeros(2,length(sigma2s));ref_Es=zeros(1,length(sigma2s));
srt_E0=zeros(2,length(sigma2s));ref_E0=zeros(2,length(sigma2s));

repeat_max=1;
for i=(1:length(sigma2s))
    sigma2=sigma2s(i)
    C_qos=abs(C_qos0/log(10^-14/sigma2)*log(2));
    if sigma2==10^-14
        C_qos=C_qos0;
    end
    h0=ones(J_range+1,J_range,T_range);
    [srt_E0(1,i),srt_E0(2,i),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
    [ref_E0(1,i),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
    srt_Es_partial=zeros(2,repeat_max);ref_Es_partial=zeros(1,repeat_max);
    for ii=(1:repeat_max)
        
        h0=normrnd(0,1,[J_range+1,J_range,T_range*h0_multi]);
        [srt_Es_partial(1,ii),srt_Es_partial(2,ii),~]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);
        [ref_Es_partial(1,ii),~]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0);    
    end
    srt_Es(:,i)=mean(srt_Es_partial,2);ref_Es(:,i)=mean(ref_Es_partial,2);
    srt_E0(:,i)=srt_E0(:,i)/J_range/C_qos(1);ref_E0(:,i)=ref_E0(:,i)/J_range/C_qos(1);
    srt_Es(:,i)=srt_Es(:,i)/J_range/C_qos(1);ref_Es(:,i)=ref_Es(:,i)/J_range/C_qos(1); 
end

save('.\experiments\sigma2s','sigma2s','srt_Es','ref_Es','srt_E0','ref_E0') 
plot(sigma2s,srt_E0(1,:),'r--',sigma2s,srt_E0(2,:),'r--o',sigma2s,ref_E0(1,:),'g--',sigma2s,srt_Es(1,:),'r-',sigma2s,srt_Es(2,:),'r-o',sigma2s,ref_Es(1,:),'g-');
xlabel('sigma^2')
ylabel('Energy Consumption per User per Gbit (J)')
legend('h_0 constant, cellular-only','h_0 constant, D2D underlaid','h_0 constant, reference','h_0 gaussian, cellular-only','h_0 gaussian, D2D underlaid','h_0 gaussian, reference')

