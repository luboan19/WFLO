clear;yalmip('clear');
address=[];
for x=0:200:1800
    for y=0:200:1800
         address=[address;[x y]];
    end
end
% wind=[12,0,1];
wind=[];
for n=0:11
    wind=[wind;[12,2*pi/12*n,1/12]];
end
DD=Wake_effect_v1(address,wind);
WW=DD.^2;
I=size(address,1); %I为将要选取的地址的离散集合,格式为address(x1 y1;x2 y2...）
D=size(wind,1);%D为风的离散集合,格式为wind(v1 ang1 p1;v2 ang2 p2;.....），角度为与正北方向按顺时针计算的角度大小
y=binvar(size(address,1),1,'full');
loss_m=sdpvar(size(address,1),size(wind,1),'full');
loss_sqrt=sdpvar(size(address,1),size(wind,1),'full');
v_rate=sdpvar(size(address,1),size(wind,1),'full');
v_fang=sdpvar(size(address,1),size(wind,1),'full');
v_san=sdpvar(size(address,1),size(wind,1),'full');
p=sdpvar(size(address,1),size(wind,1),'full');
Constraints=[sum(y)==30,p>=0,v_rate>=0,v_fang>=0,v_san>=0,loss_m>=0,loss_sqrt>=0];  %约束，总机数
for dd=1:D
    Constraints=[Constraints,p(:,dd)<=0.3*wind(dd,3)*(wind(dd,1)^3)*y];
end
tic
for dd=1:D
    for jj= 1:I
    Constraints=[Constraints,loss_m(jj,dd) >=y'*WW(:,jj,dd)];    
    Constraints=[Constraints,loss_sqrt(jj,dd)*loss_sqrt(jj,dd) >=loss_m(jj,dd)];    
    Constraints=[Constraints,v_rate(jj,dd) <=(1-loss_sqrt(jj,dd))]; 
    Constraints=[Constraints,v_fang(jj,dd)<=v_rate(jj,dd)*v_rate(jj,dd)];   
    Constraints=[Constraints,v_san(jj,dd)<=v_rate(jj,dd)*v_fang(jj,dd)];  
    Constraints=[Constraints,p(jj,dd)<=0.3*wind(dd,3)*(wind(dd,1)^3)*v_san(jj,dd)];
    end
end
toc
Objective=sum(p,'all');
ops = sdpsettings('solver', 'gurobi','verbose', 2,'gurobi.Timelimit',600,'usex0',1);%设置求解器
% read= readmatrix('quarcombine.txt');  %赋初值
% assign(y,read);
solution=optimize(Constraints,-Objective,ops)
aaa=value(y);
for n=1:length(aaa)
    if((0.999<=aaa(n)&&aaa(n)<=1.001))
        plot(address(n,1),address(n,2),'pentagram');
        hold on;
    end
end
title('风机微观选址布局方案');xlabel('横轴位置/m');ylabel('纵轴位置/m');
writematrix(aaa,'quarcombine.txt')
loss=0;
turbine_number=0;
power=0;
for dd=1:D
    for jj= 1:I
        for ii=1:I
            loss=loss+(DD(ii,jj,dd)^2)*aaa(ii);
        end
        turbine_number=turbine_number+aaa(jj)*(1-sqrt(loss))^3;  %用v^3计算功率
        loss=0;
    end
    power=power+0.3*wind(dd,3)*(wind(dd,1)^3)*turbine_number;
    turbine_number=0;
end
% for dd=1:D
%     for jj= 1:I
%         for ii=1:I
%             loss=loss+(1-vv(ii,jj,dd)^3)*aaa(ii)*aaa(jj);
%         end
%         turbine_number=turbine_number+aaa(jj)-loss;
%         loss=0;
%     end
%     power=power+0.3*wind(dd,3)*(wind(dd,1)^3)*turbine_number;
%     turbine_number=0;
% end
power
p_obj=value(Objective)


