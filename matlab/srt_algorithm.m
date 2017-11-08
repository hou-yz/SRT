function [E_stage_2,E_stage_3,C_unsatisfied]=srt_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et)
% service begin time -- bt
% service end time   -- et
P=reshape([P_i_max;ones(J_range,1)*P_j_max]*ones(1,J_range*T_range),J_range+1,J_range,T_range);
P_s=linspace(0,P_i_max,101);
P_i_tmp=ones(J_range+1,J_range,T_range,1).*reshape(P_s,[1,1,1,101]);
%fast fading, size(h0,3)>T_range
T_range_h0=size(h0,3);
t_step=T_range_h0/T_range;
r=zeros(J_range+1,J_range,T_range,101);
F_x= @(x,z) exp(-z.*x)./x;
if t_step>1
    for t=(1:T_range)
        t_h0=(t-1)*t_step;
        r(:,:,t,:)=sum(BW*log2(1+P_i_tmp(:,:,t,:).*(beta(:,:,t).*ones(1,1,t_step,101)).*(abs(h0(:,:,t_h0+1:t_h0+t_step)).^2.*ones(1,1,1,101))/sigma2),3);
    end
    r=r/t_step;
else
    gamma=beta.*ones(1,1,1,101).*P_i_tmp./sigma2;
    F_int=integral(@(x) F_x(x,1./gamma),1,inf,'ArrayValued',true);
    r=BW*log2(exp(1)).*exp(1./gamma).*F_int;
    r(isnan(r))=0;r(r==inf)=0;
end

r_max=[r(1,:,:,101);r(2:J_range+1,:,:,11)];

S=[];%[i,j,t]->[0,j,t,0,0,0]
%以下为算法第一部分
% eta = 传输时间 / 时隙长度
eta=zeros(J_range+1,J_range,T_range);
C_tmp=zeros(J_range,T_range+1);
C_plausible=C_qos;
%C_tmp(j,t) means have that much data at the time slot before t

for j=(1:J_range) %遍历用户j
    [~,beta_j_sorted_t_index]=sort(r(1,j,bt(j):et(j),101),'descend');
    ijt_t_pointer=1;
    while C_tmp(j,T_range+1)<C_plausible(j) %若不满足qos条件，增加一条边
        if isempty(beta_j_sorted_t_index)
            C_plausible(j)=C_tmp(j,T_range+1);
            break
        end
        t0=beta_j_sorted_t_index(ijt_t_pointer)+bt(j)-1;
        if r(0+1,j,t0,101)<=0
            %fprintf('srt:用户%d未能满足QoS限制',j)
            C_plausible(j)=C_tmp(j,T_range+1);
            break
        end
        if C_tmp(j,T_range+1)+r(0+1,j,t0,101)*delta_t<C_qos(j)%未溢出
            eta(0+1,j,t0)=1;%index 0 for BS-i
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r(0+1,j,t0,101)*delta_t;
        else%溢出
            r0=(C_qos(j)-C_tmp(j,T_range+1))/delta_t;
%             [~,P0_index]=min(abs(r(0+1,j,t0,:)-r0));
%             P0=P_s(P0_index);
%             eta(0+1,j,t0)=P0/P_i_max;
            eta(0+1,j,t0)=r0/r(0+1,j,t0,101);
            C_tmp(j,t0:T_range+1)=C_plausible(j);
        end
        
        S=[S;[0,j,t0,0,0,0]];
        ijt_t_pointer=ijt_t_pointer+1;
        if ijt_t_pointer>et(j)-bt(j)+1
            %fprintf('srt:用户%d未能满足QoS限制',j)
            C_plausible(j)=C_tmp(j,T_range+1);
        end
    end
    
end
%第一部分结束
E_stage_1=sum(sum(sum(P.*eta*delta_t)));
%fprintf('第一部分结束，系统总功耗：%d\n',E_stage_1)
delta_1=reshape(eta(0+1,:,:),J_range,T_range);

%以下为算法第二部分
t_free=ones(1,T_range);%time slots that are availible 初始全部时隙为空闲
t_OOC=zeros(1,T_range);%out of capacity 初始无时隙超出服务能力
for t=(1:T_range)
    t_free(t)=(sum(eta(0+1,:,t)>0)<N);
    t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
end
while sum(t_OOC)~=0 %超出基站能力范围
    %更新
    for t=(1:T_range)
        t_free(t)=(sum(eta(0+1,:,t)>0)<N);
        t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
    end
    r_i_free=reshape(r(1,:,:).*(1-eta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free.*(ones(J_range,1)*(1:T_range)>=reshape(bt,[J_range,1])*ones(1,T_range)).*(ones(J_range,1)*(1:T_range)<=reshape(et,[J_range,1])*ones(1,T_range)));% 空闲时隙中 可用 速率
    r_i_OOC=reshape(r(1,:,:),J_range,T_range).*(ones(J_range,1)*t_OOC);r_i_OOC(eta(0+1,:,:)==0)=0;r_i_OOC(r_i_OOC==0)=NaN;% 超出能力时隙对应速率
    delta_C=zeros(1,J_range);% 记录对用户j，容量变化（对系统功率影响）最小解
    t0_index=zeros(1,J_range);% 每个用户j，容量变化（对系统功率影响）最小解对应 原 时隙t0
    for j=(1:J_range)
        [delta_C(j),t0_index(j)]=min(r_i_OOC(j,:)*delta_t);
        delta_C(j)=delta_C(j)-max(r_i_free(j,:))*delta_t;
    end
    [~,j]=min(delta_C); %去掉效率（beta）较差的边
    t0=t0_index(j);
    
    [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
    S(S_index,:)=[];
    C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r(0+1,j,t0)*eta(0+1,j,t0)*delta_t;
    eta(0+1,j,t0)=0;
    %fprintf('去掉一条边：[0 %d %d]\n',j,t0)
    %更新
    for t=(1:T_range)
        t_free(t)=(sum(eta(0+1,:,t)>0)<N);
        t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
    end
    r_i_free=reshape(r(1,:,:).*(1-eta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free.*(ones(J_range,1)*(1:T_range)>=reshape(bt,[J_range,1])*ones(1,T_range)).*(ones(J_range,1)*(1:T_range)<=reshape(et,[J_range,1])*ones(1,T_range)));% 空闲时隙中 可用 速率
    
    while C_tmp(j,T_range+1)<C_plausible(j)%去掉该边后，寻找一组次优解满足QoS限制
        if sum(r_i_free(j,:))<=0
            %fprintf('srt:用户%d未能满足QoS限制',j)
            C_plausible(j)=C_tmp(j,T_range+1);
            break;
        end
        [r_ijt0,t0]=max(r_i_free(j,:));
        if C_tmp(j,T_range+1)+r_ijt0*delta_t<C_plausible(j)%未溢出
            eta(0+1,j,t0)=1;% 和之前delta取值无关，现在全都更新为1
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r_ijt0*delta_t*eta(0+1,j,t0);
        else%溢出
            eta(0+1,j,t0)=eta(0+1,j,t0)+(C_plausible(j)-C_tmp(j,T_range+1))/(r_ijt0*delta_t);
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)+r_ijt0*delta_t*eta(0+1,j,t0);
        end
        
        S=[S;[0,j,t0,0,0,0]];
        %fprintf('增加一条边：[0 %d %d]\n',j,t0)
        %更新
        for t=(1:T_range)
            t_free(t)=(sum(eta(0+1,:,t)>0)<N);
            t_OOC(t)=(sum(eta(0+1,:,t)>0)>N);
        end
        r_i_free=reshape(r(1,:,:).*(1-eta(0+1,:,:)),J_range,T_range).*(ones(J_range,1)*t_free.*(ones(J_range,1)*(1:T_range)>=reshape(bt,[J_range,1])*ones(1,T_range)).*(ones(J_range,1)*(1:T_range)<=reshape(et,[J_range,1])*ones(1,T_range)));% 空闲时隙中 可用 速率
    
    end    
end 
%第二部分结束
E_stage_2=sum(sum(sum(P.*eta*delta_t)));
%fprintf('第二部分结束，系统总能耗：%d\n',E_stage_2)
%reshape(delta(0+1,:,:),J_range,T_range)
delta_2=reshape(eta(0+1,:,:),J_range,T_range);
C_unsatisfied=C_qos-C_plausible;

%第三部分
r0_index=find(eta(0+1,:,:)>0);
[~,jj,tt]=ind2sub(size(eta(0+1,:,:)),r0_index);
for j=(1:J_range)%for all user j
    R=[];
    subs_found=0;
    r0_j=r(0+1,j,tt(jj==j),101);
    if isempty(r0_j)
        continue
    end
    [r0_max,~]=max(r0_j);
    for j_prime_R=(1:J_range)
        
        t2s=(bt(j):et(j));
        t2s=t2s(r(j_prime_R+1,j,t2s,11)>r0_max);
        for t2=t2s
            %disp([j_prime_R j sum(delta(:,j_prime_R,t2)>0)+sum(delta(j_prime_R+1,:,t2)>0)-(delta(j_prime_R+1,j,t2)>0) sum(delta(:,j,t2)>0)-(delta(j_prime_R+1,j,t2)>0)+sum(delta(j+1,:,t2)>0) t2])
            if j_prime_R~=j && sum(eta(:,j_prime_R,t2)>0)+sum(eta(j_prime_R+1,:,t2)>0)-2*(eta(j_prime_R+1,j,t2)>0)==0 && sum(eta(:,j,t2)>0)-2*(eta(j_prime_R+1,j,t2)>0)+sum(eta(j+1,:,t2)>0)==0 && sum(sum(eta(:,:,t2)>0))-(eta(j_prime_R+1,j,t2)>0)<N && eta(j_prime_R+1,j,t2)<1% j' and j and sys free
                
                t1s=(bt(j):et(j));
                t1s=t1s(r(0+1,j_prime_R,t1s,101)>r0_max);
                for t1=t1s
                    if t1~=t2 && sum(eta(:,j_prime_R,t1)>0)-2*(eta(0+1,j_prime_R,t1)>0)+sum(eta(j_prime_R+1,:,t1)>0)==0 && sum(sum(eta(:,:,t1)>0))-(eta(0+1,j_prime_R,t1)>0)<N && eta(0+1,j_prime_R,t1)<1 % j' free and sys availible
                        r1_max=r(0+1,j_prime_R,t1,101);%r1
                        if ((t1>t2 && C_tmp(j_prime_R,t2)>=r0_max*delta_t) || (t1<=t2 && C_tmp(j_prime_R,t2)+r1_max*delta_t>=r0_max*delta_t))
                              R=[R;[0,j_prime_R,t1,j_prime_R,j,t2]];
                        end
                    end
                end
            end
        end
        
    end
    delta_R=zeros(size(R,1),1);
    for index=(1:size(R,1))
        delta_R(index,:)=10/r_max(R(index,1)+1,R(index,2),R(index,3))+1/r_max(R(index,4)+1,R(index,5),R(index,6));
    end
    [~,index_Rs]=sort(delta_R);
    t0s=find(eta(0+1,j,:)>0);
    r0s=reshape((r_max(0+1,j,t0s)),[length(t0s),1]);
    [~,index_0s]=sort(r0s);
    i_R_pointer=1;i_S_pointer=1;
    C_r_1=0;C_r_2=0;
    while i_R_pointer<=size(R,1) && i_S_pointer<=length(r0s)
        index_R=index_Rs(i_R_pointer);
        j_prime=R(index_R,2);t1=R(index_R,3);t2=R(index_R,6);
        r1=r(0+1,j_prime,t1,101);r2=r(j_prime+1,j,t2,11);
        index_0=index_0s(i_S_pointer);
        t0=t0s(index_0);
        r0=r0s(index_0);
        eta0=eta(0+1,j,t0);
        if C_r_1 + C_r_2 == 0
            eta1=r0*eta(0+1,j,t0)/r1;
            eta2=r0*eta(0+1,j,t0)/r2;
        else
            eta1=C_r_1/r1/delta_t;
            eta2=C_r_2/r2/delta_t;
        end
        if sum(eta(:,j_prime,t1)>0)-(eta(0+1,j_prime,t1)>0)+sum(eta(j_prime+1,:,t1)>0)==0 && sum(eta(:,j_prime,t2)>0)+sum(eta(j_prime+1,:,t2)>0)-(eta(j_prime+1,j,t2)>0)==0 && sum(eta(j+1,:,t2)>0)+sum(eta(:,j,t2)>0)-(eta(j_prime+1,j,t2)>0)==0 && sum(sum(eta(:,:,t2)>0))-(eta(j_prime+1,j,t2)>0)<N && sum(sum(eta(:,:,t1)>0))-(eta(0+1,j_prime,t1)>0)<N
            %eta(0+1,j_prime,t1)+eta1<1 && eta(j_prime+1,j,t2)+eta2<1 && 
            i_R_p_last=0;i_S_p_last=0;
        else
            if ~i_R_p_last && ~i_S_p_last
                %初次出现问题，记录
                i_R_p_last=i_R_pointer;i_S_p_last=i_S_pointer;
            else
                %已经出现过问题，循环遍历直到找到可行解
                i_R_pointer=i_R_pointer+1;
                if i_R_pointer==size(R,1) 
                    i_S_pointer=i_S_pointer+1;
                    i_R_pointer=i_R_p_last;
                end
            end
            continue
        end
        if eta0*10-eta1*10-eta2>0
            if C_r_1 + C_r_2 == 0
                %去除edge0边.
                eta(0+1,j,t0)=0;
                [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
                if S_index~=0
                    S(S_index,:)=[];
                    if ~i_S_p_last
                        i_S_pointer=i_S_pointer+1;
                    else 
                        i_S_pointer=i_S_p_last+1;
                    end
                else
                    break
                end
                C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r0*delta_t;
                %fprintf('去掉一条边：[0 %d %d]\n',j,t0)
            end
            %可以满足的eta1 eta2
            eta1_s=min((1-eta(0+1,j_prime,t1)),eta1);eta2_s=min((1-eta(j_prime+1,j,t2)),eta2);
            %尚未满足的eta1 eta2
            C_r_1=(eta1-eta1_s)*r1*delta_t;C_r_2=(eta2-eta2_s)*r2*delta_t;
            
            %加入edge2边
            eta(j_prime+1,j,t2)=eta(j_prime+1,j,t2)+eta2_s;
            C_tmp(j,t2:T_range+1)=C_tmp(j,t2:T_range+1)+r2*eta2_s*delta_t;
            C_tmp(j_prime,t2:T_range+1)=C_tmp(j_prime,t2:T_range+1)-r2*eta2_s*delta_t;%*(1-gamma(j_prime,j))
            %加入edge1边
            eta(0+1,j_prime,t1)=eta(0+1,j_prime,t1)+eta1_s;
            S=[S;[0,j_prime,t1,j_prime,j,t2]];
            C_tmp(j_prime,t1:T_range+1)=C_tmp(j_prime,t1:T_range+1)+r1*eta1_s*delta_t;%*(1-gamma(j_prime,j))
            %fprintf('增加两条边：[0 %d %d] [%d %d %d]\n',j_prime,t1,j_prime,j,t2)
                
            if C_r_1 + C_r_2
                if ~i_R_p_last
                    i_R_pointer=i_R_pointer+1;
                else 
                    i_R_pointer=i_R_p_last+1;
                end
            end
            
            
        else
            break
        end
    
    end
    
    %{
    if i_R_pointer<size(R,1) && i_S_pointer<length(r0s)
    while ~isempty(R)
        delta_E=zeros((length(t0s)-i_S_pointer)*(size(R,1)-i_R_pointer),3);edges=zeros((length(t0s)-i_S_pointer)*(size(R,1)-i_R_pointer),9);
        for i_0=(i_S_pointer:length(t0s))
            index_0=index_0s(i_0);
            t0=t0s(index_0);
            r0=r_max(0+1,j,t0);
            edge0=[0,j,t0];
%             r0=r(0+1,j,t0,round(eta(0+1,j,t0)*100+1));
%             P0=P_i_max*eta(0+1,j,t0);
            P0=P_i_max*eta(0+1,j,t0);

            for i_12=(i_R_pointer:size(R,1))
                index_R=index_Rs(i_12);
                j_prime=R(index_R,2);t1=R(index_R,3);t2=R(index_R,6);
                r1=r(0+1,j_prime,t1,101);r2=r(j_prime+1,j,t2,11);
%                 [~,P1_index]=min(abs(r(0+1,j_prime,t1,:)-r0));
%                 P1=P_s(P1_index);
%                 eta1=P1/P_i_max;
%                 [~,P2_index]=min(abs(r(j_prime+1,j,t2,:)-r0));
%                 P2=P_s(P2_index);
%                 eta2=P2/P_i_max;
                eta1=r0*eta(0+1,j,t0)/r(0+1,j_prime,t1,101);
                eta2=r0*eta(0+1,j,t0)/r(j_prime+1,j,t2,11);
                P1=P_i_max*eta1;P2=P_j_max*eta2;
                ii=(i_0-1)*size(R,1)+i_12;
                if eta(0+1,j_prime,t1)+eta1<1 && eta(j_prime+1,j,t2)+eta2<1 && sum(eta(:,j_prime,t1)>0)-(eta(0+1,j_prime,t1)>0)+sum(eta(j_prime+1,:,t1)>0)==0 && sum(eta(:,j_prime,t2)>0)+sum(eta(j_prime+1,:,t2)>0)-(eta(j_prime+1,j,t2)>0)==0 && sum(eta(j+1,:,t2)>0)+sum(eta(:,j,t2)>0)-(eta(j_prime+1,j,t2)>0)==0 && sum(sum(eta(:,:,t2)>0))-(eta(j_prime+1,j,t2)>0)<N && sum(sum(eta(:,:,t1)>0))-(eta(0+1,j_prime,t1)>0)<N
                    delta_E(ii,:)=[eta1,eta2,(P0-(P1+P2))*delta_t];
                else
                    delta_E(ii,:)=[eta1,eta2,-inf];
                end
                edges(ii,:)=[edge0,edge1,edge2];
            end
        end
    
        [m,ii]=max(delta_E(:,3));%得到最优的替换方式 delta_P(ii,:)=[P1,P2,P_i_max-(P1+P2)];
        if m>0
            %替换
            edge0=edges(ii,1:3);edge1=edges(ii,4:6);edge2=edges(ii,7:9);
            j_prime=edge1(2);
            t0=edge0(3);t1=edge1(3);t2=edge2(3);
            r0=r(0+1,j,t0,round(eta(0+1,j,t0)*100+1));
            eta1=delta_E(ii,1);eta2=delta_E(ii,2);
            %去除edge0边.
            eta(0+1,j,t0)=0;
            [~,S_index]=ismember([0,j,t0,0,0,0],S,'rows');
            if S_index~=0
                S(S_index,:)=[];
            else
                break
            end
            C_tmp(j,t0:T_range+1)=C_tmp(j,t0:T_range+1)-r0*delta_t;
            %fprintf('去掉一条边：[0 %d %d]\n',j,t0)
            %加入edge2边
            eta(j_prime+1,j,t2)=eta(j_prime+1,j,t2)+eta2;
            C_tmp(j,t2:T_range+1)=C_tmp(j,t2:T_range+1)+r0*delta_t;
            C_tmp(j_prime,t2:T_range+1)=C_tmp(j_prime,t2:T_range+1)-r0*delta_t;%*(1-gamma(j_prime,j))
            %加入edge1边
            eta(0+1,j_prime,t1)=eta(0+1,j_prime,t1)+eta1;
            S=[S;[0,j_prime,t1,j_prime,j,t2]];
            C_tmp(j_prime,t1:T_range+1)=C_tmp(j_prime,t1:T_range+1)+r0*delta_t;%*(1-gamma(j_prime,j))
            %fprintf('增加两条边：[0 %d %d] [%d %d %d]\n',j_prime,t1,j_prime,j,t2)
            
            %从R中去除
            %[~,R_index]=ismember([0,j_prime,t1,j_prime,j,t2],R,'rows');
%             tmp_max=size(R,1);tmp=1;
%             while tmp<tmp_max
%                 if R(tmp,6)==t2 || ( R(tmp,2)==j_prime && R(tmp,3)==t1 )
%                     R(tmp,:)=[];
%                     tmp=tmp-1;
%                     tmp_max=tmp_max-1;
%                 end
%                 tmp=tmp+1;
%             end
        else
            break
        end
    
    end
    end
    %}
end
E_stage_3=sum(sum(sum(P.*eta*delta_t)));
%fprintf('第三部分结束，系统总能耗：%d\n',E_stage_3)

end    