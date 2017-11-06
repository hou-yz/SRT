function [E_stage_1,C_unsatisfied]=ref_algorithm(beta,J_range,T_range,P_i_max,P_j_max,C_qos,N,delta_t,BW,sigma2,h0,bt,et)
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


S=[];%[i,j,t]->[0,j,t,0,0,0]
C_max=C_qos;
%以下为算法第一部分
eta=zeros(J_range+1,J_range,T_range);
C_tmp=zeros(J_range,T_range+1);
%C_tmp(j,t) means have that much data at the time slot before t
satisfied=zeros(1,J_range);
j=1;
for t0=(1:T_range) %遍历时隙
    j_add=0;
    if prod(satisfied)
        break
    end
    while sum(C_tmp(:,T_range+1)<C_qos(:)) && sum(eta(1,:,t0)>0)<N && j_add<J_range %若不满足qos条件，增加一条边  
        if t0<bt(j) || t0>et(j) || satisfied(j)
            j=j+1;
            j_add=j_add+1;
            
            if j>J_range
                j=j-J_range;
            end
            continue
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
            satisfied(j)=1;
            C_tmp(j,t0:T_range+1)=C_max(j);
        end
        S=[S;[0,j,t0,0,0,0]];
        j=j+1;
        j_add=j_add+1;
        if j>J_range
            j=j-J_range;
        end
        if prod(satisfied)
            break
        end
    end
end
if ~prod(satisfied)
    for j=(1:J_range)
        if ~satisfied(j)
        %fprintf('ref:用户%d未能满足QoS限制',j)
            C_max(j)=C_tmp(j,T_range+1);
        end
    end
end
%第一部分结束
E_stage_1=sum(sum(eta(1,:,:)))*P_i_max*delta_t;
%fprintf('第一部分结束，系统总功耗：%d\n',E_stage_1)
delta_1=reshape(eta(0+1,:,:),J_range,T_range);
C_unsatisfied=C_qos-C_max;

end

