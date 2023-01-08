  function DD = Wake_effect(address,wind)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �ڴ�ģ������ʱ�������ص������Ӱ�죬ʹ�÷Ǹ�˹jessenģ��
I=size(address,1); %IΪ��Ҫѡȡ�ĵ�ַ����ɢ����,��ʽΪ(x1 y1;x2 y2...��
D=size(wind,1);%DΪ�����ɢ����,��ʽΪ(v1 ang1;v2 ang2;.....�����Ƕ�Ϊ����������˳ʱ�����ĽǶȴ�С
alpha=0.1;
rd=20;
DD=zeros(I,I,D);
Ct=0.88;
    for k=1:D  %ÿ�����򡢷���
        for j=1:I    %ÿ����Ӱ���j
            for i=1:I   %��ÿ��j�����˵��ͬ��i��������Ӱ��
                d=(address(i,2)-address(j,2))*cos(wind(k,2))+(address(i,1)-address(j,1))*sin(wind(k,2));%dΪ��������ڷ����ϵľ���
                d_v=abs((address(i,2)-address(j,2))*sin(wind(k,2))-(address(i,1)-address(j,1))*cos(wind(k,2)));%lΪ��������ڷ����ϵĴ�ֱ����
                if(d<=0)  %dС�ڵ���0��β��������Ӱ��
                    DD(i,j,k)=0;
                else     %d>0
                    if(wind(k,2)==pi/2||wind(k,2)==3*pi/2)   %�Ƕ�Ϊ��������
                        y_max=address(i,2)-d*cos(wind(k,2))+(rd+d*alpha+rd)*abs(sin(wind(k,2)));
                        y_min=address(i,2)-d*cos(wind(k,2))-(rd+d*alpha+rd)*abs(sin(wind(k,2)));
                        if(y_min<=address(j,2)&&address(j,2)<=y_max)  %��β��ЧӦӰ��
                            if(d_v<=d*alpha)
                            DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2);
                            elseif(d*alpha<d_v&&d_v<=rd)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+(pi-acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v)))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v
                               ita=A_shadow/(pi*rd^2)
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita; 
                            elseif(d_v>rd)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v
                               ita=A_shadow/(pi*rd^2)
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita;   
                            end
                        else
                            DD(i,j,k)=0;
                        end
                    else  %�Ƕȷ���������
                        x_max=address(i,1)-d*sin(wind(k,2))+(rd+d*alpha+rd)*abs(cos(wind(k,2)));
                        x_min=address(i,1)-d*sin(wind(k,2))-(rd+d*alpha+rd)*abs(cos(wind(k,2)));
                        if(x_min<=address(j,1)&&address(j,1)<=x_max)  %��β��ЧӦӰ��
                            if(d_v<=d*alpha)
                            DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2);   %�Ľ�jessenβ��ģ��
                            elseif (d*alpha<d_v&&d_v<=rd) 
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+(pi-acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*(rd)*d_v)))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v
                               ita=A_shadow/(pi*rd^2);
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita;  %�Ľ�jessenβ��ģ��
                            elseif(d_v>rd)  
                               A_shadow=acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v))*(d*alpha+rd)^2+acos((-(d*alpha+rd)^2+d_v^2+rd^2)/(2*rd*d_v))*(rd)^2-sin(acos(((d*alpha+rd)^2+d_v^2-rd^2)/(2*(d*alpha+rd)*d_v)))*(d*alpha+rd)*d_v
                               ita=A_shadow/(pi*rd^2);
                               DD(i,j,k)=(1-sqrt(1-Ct))*((rd/(rd+alpha*d))^2)*ita;  %�Ľ�jessenβ��ģ�� 
                            end                                
                        else
                            DD(i,j,k)=0;
                        end
                    end
                end
            end
        end
    end
end

