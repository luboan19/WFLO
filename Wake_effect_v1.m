   function DD = Wake_effect_v1(address,wind)
%   β����ʧ���ɺ��� �����ʽΪ(address,wind)
%   ���о����ʽaddress(x1 y1;x2 y2...) wind(v1 ang1 p1;v2 ang2 p2;.....) angΪ����������˳ʱ����ת����ĽǶ� 
%   �ڴ�ģ���п����ص������Ӱ�죬ʹ�÷Ǹ�˹�ֲ���jessenģ��
I=size(address,1); %IΪ��Ҫѡȡ�ĵ�ַ����ɢ����,��ʽΪ(x1 y1;x2 y2...��
D=size(wind,1);%DΪ�����ɢ����,��ʽΪ(v1 ang1 p1 ;v2 ang2 p2;.....�����Ƕ�Ϊ����������˳ʱ�����ĽǶȴ�С
alpha=0.1;
rd=20;
DD=zeros(I,I,D);
Ct=0.88;
    for k=1:D  %ÿ�����򡢷���
        for j=1:I    %ÿ����Ӱ���j
            for i=1:I   %��ÿ��j�����˵��ͬ��i��������Ӱ��
                d=(address(i,2)-address(j,2))*cos(wind(k,2))+(address(i,1)-address(j,1))*sin(wind(k,2));%dΪ��������ڷ����ϵľ���
                d_v=abs((address(i,2)-address(j,2))*sin(wind(k,2))-(address(i,1)-address(j,1))*cos(wind(k,2)));%lΪ��������������ֱ�����ϵľ���
                if(d<=0)  %dС�ڵ���0��β��������Ӱ��
                    DD(i,j,k)=0;
                else     %d>0
                            if(d_v<=d*alpha)
                            DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2);
                            elseif(d*alpha<d_v&&d_v<=rd+d*alpha)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+(pi-acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v)))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v;
                               ita=A_shadow/(pi*rd^2);
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita; 
                            elseif(d_v>rd+d*alpha&&d_v<=rd+rd+d*alpha)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v;
                               ita=A_shadow/(pi*rd^2);
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita;   
                            elseif(d_v>rd+rd+d*alpha)  
                            DD(i,j,k)=0;
                            end
                end
            end
        end
    end
end

