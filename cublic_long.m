    %-------------------------------------------------------------------
%--------     cublic function  nagtive defination  long   ----------------
%--------     Fixed Tangent -> Dynamic Tangent   ----------------
clc
clear
close all
tic

%% 系统矩阵设置
%系统维数
n = 2;
%可调参数
namuda = 1;
%给定时延下界
h1 = 1; 

%二分查找范围
ACC1 = h1;
ACC2 = 25; 
acc = 0.001;

while true
h2 = (ACC1 + ACC2) / 2;
h12 = h2 - h1; 
y0=(1-namuda) * h1 + namuda * h2;

%% LMI变量定义
num=12;
for i = 1 : num
    e = [zeros(n,(i-1) * n), eye(n), zeros(n,(num-i) * n)]; % 定义ei
    eval(['e',num2str(i),'=[e]']); %计算 MATLAB 表达式，即计算e
end
A = [ 0, 1;
     -10, -1];
Ad = [0, 0.1;
      0.1, 0.2];
es=A*e1+Ad*e3; 
setlmis([])
[P,  P_n,   sP] = lmivar(1,[7*n,1]);
[Q1, Q1_n, sQ1] = lmivar(1,[2*n,1]);
[Q2, Q2_n, sQ2] = lmivar(1,[6*n,1]);
[R1, R1_n, sR1] = lmivar(1,[n,1]);
[R2, R2_n, sR2] = lmivar(1,[n,1]);
[N1, N1_n, sN1] = lmivar(1,[3*n,1]);
[N5, N5_n, sN5] = lmivar(1,[3*n,1]);
[N2, N2_n, sN2] = lmivar(2,[3*n,3*n]);
[N4, N4_n, sN4] = lmivar(2,[3*n,3*n]);
[N3, N3_n, sN3] = lmivar(2,[3*n,5*n]);
hatR1 = lmivar(3,[sR1,        zeros(n,n),  zeros(n,n);
                zeros(n,n), sR1,         zeros(n,n);
                zeros(n,n), zeros(n,n),       sR1]);
hatR2 = lmivar(3,[sR2,        zeros(n,n),  zeros(n,n);
                zeros(n,n), sR2,         zeros(n,n);
                zeros(n,n), zeros(n,n),  sR2]);
L1=[e1;e2;e4;h1*e5;-h1*e6+h2*e7;h1*h1*e8;h1*h1*e9+h2*h2*e10-h1*h2*e6];
L2=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);e6-e7;zeros(n,num*n);-2*h1*e9-2*h2*e10+(h1+h2)*e6];
L3=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);e9+e10-e6];
L4=[es;e11;e12;e1-e2;e2-e4;h1*(e1-e5);h12*e2+h1*e6-h2*e7];
L5=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);-e6+e7];
L6=[e2;e11;e1;e2;e4;zeros(n,num*n)];
L7=[e4;e12;e1;e2;e4;-h1*e6+h2*e7];
L8=[-h1*e6+h2*e7;e2-e4;h12*e1;h12*e2;h12*e4;h1*h1*e9+h2*h2*e10-h1*h2*e6];
L9=[e6-e7;zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);-2*h1*e9-2*h2*e10+(h1+h2)*e6];
L10=[zeros(n,num*n);zeros(n,num*n);es;e11;e12;e2];
L11=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);e9+e10-e6];
L12=[zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);zeros(n,num*n);-e6+e7];
gamma0=[e1;e2;e3;e11;e12];
gamma1=[e1-e2; e1+e2-2*e5; e1-e2+6*e5-12*e8];
gamma2=[e2-e3; e2+e3-2*e6; e2-e3+6*e6-12*e9];
gamma3=[e3-e4; e3+e4-2*e7; e3-e4+6*e7-12*e10];
Coe=[eye(n),      zeros(n),         zeros(n);
    zeros(n),    sqrt(3)*eye(n),   zeros(n);
    zeros(n),    zeros(n),         sqrt(5)*eye(n)];
%%
i=1;
lmiterm([i,1,1,P],L1'+h1*L2'+h1*h1*L3',L4+h1*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h1*L12',L7-h1*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h1*L9'+h1*h1*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-2*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N1],-gamma2',gamma2,'s');
lmiterm([i,1,1,N2],-gamma2',gamma3,'s');
lmiterm([i,1,1,N3],-gamma2',gamma0,'s');

lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);

%%
i=2;
lmiterm([i,1,1,P],L1'+h2*L2'+h2*h2*L3',L4+h2*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h2*L12',L7-h2*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h2*L9'+h2*h2*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-2*gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N4],-gamma3',gamma2,'s');
lmiterm([i,1,1,N5],-gamma3',gamma3,'s');
lmiterm([i,1,1,N3],gamma3',gamma0,'s');

lmiterm([i,1,2,-N1],gamma2',1);
lmiterm([i,1,2,-N2],gamma3',1);
lmiterm([i,1,2,-N3],gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);
%% 
i=3;
lmiterm([i,1,1,P],L1'+h1*L2'+h1*h1*L3',L4+h1*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h1*L12',L7-h1*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h1*L9'+h1*h1*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-2*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N1],-gamma2',gamma2,'s');
lmiterm([i,1,1,N2],-gamma2',gamma3,'s');
lmiterm([i,1,1,N3],-gamma2',gamma0,'s');
lmiterm([i,1,1,P],-3*h12*h12*h2*L3',L5,'s');
lmiterm([i,1,1,P],-h12*h12*L2',L5,'s');
lmiterm([i,1,1,P],-h12*h12*L3',L4,'s');
lmiterm([i,1,1,Q2],-h12*h12*L11',L10,'s');
lmiterm([i,1,1,Q2],h12*h12*L12',L12);

lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);

%% 
i=4;
lmiterm([i,1,1,P],L1'+h1*L2'+h1*h1*L3',L4+h1*L5,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7'+h1*L12',L7-h1*L12);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,Q2],L8'+h1*L9'+h1*h1*L11',L10,'s');
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-2*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-gamma3'*Coe,Coe*gamma3);
lmiterm([i,1,1,N1],-gamma2',gamma2,'s');
lmiterm([i,1,1,N2],-gamma2',gamma3,'s');
lmiterm([i,1,1,N3],-gamma2',gamma0,'s');
lmiterm([i,1,1,P],h12*h12*h12*L3',L5,'s');

lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);
%% 
 i=5;

lmiterm([i,1,1,P],(3*h1*y0*y0-2*y0*y0*y0)*L3',L5,'s');

lmiterm([i,1,1,P],(2*h1*y0-y0*y0)*L2',L5,'s');
lmiterm([i,1,1,P],(2*h1*y0-y0*y0)*L3',L4,'s');
lmiterm([i,1,1,Q2],(2*h1*y0-y0*y0)*L11',L10,'s');
lmiterm([i,1,1,Q2],-(2*h1*y0-y0*y0)*L12',L12);

lmiterm([i,1,1,P],h1*L2',L4,'s');
lmiterm([i,1,1,P],h1*L1',L5,'s');
lmiterm([i,1,1,Q2],h1*L9',L10,'s');
lmiterm([i,1,1,Q2],h1*L7',L12,'s');
lmiterm([i,1,1,N1],h1*(1/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],h1*(1/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],h1*(1/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],h1*(-1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],h1*(-1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],h1*(-1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,hatR2],h1*(1/h12)*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],h1*(-1/h12)*gamma3'*Coe,Coe*gamma3);
%fai_d0'
lmiterm([i,1,1,P],L1',L4,'s');
lmiterm([i,1,1,Q2],L8',L10,'s');
lmiterm([i,1,1,N1],(-h2/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],(-h2/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],(-h2/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],(h1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],(h1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],(h1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7',L7);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-(2+(h1/h12))*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-(1-(h1/h12))*gamma3'*Coe,Coe*gamma3);

%Ngamma3^T
lmiterm([i,1,2,-N4],gamma2',1);
lmiterm([i,1,2,-N5],gamma3',1);
lmiterm([i,1,2,-N3],-gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);%t>=2
%%
 i=6; 

lmiterm([i,1,1,P],(3*h2*y0*y0-2*y0*y0*y0)*L3',L5,'s');

lmiterm([i,1,1,P],(2*h2*y0-y0*y0)*L2',L5,'s');
lmiterm([i,1,1,P],(2*h2*y0-y0*y0)*L3',L4,'s');
lmiterm([i,1,1,Q2],(2*h2*y0-y0*y0)*L11',L10,'s');
lmiterm([i,1,1,Q2],-(2*h2*y0-y0*y0)*L12',L12);

lmiterm([i,1,1,P],h2*L2',L4,'s');
lmiterm([i,1,1,P],h2*L1',L5,'s');
lmiterm([i,1,1,Q2],h2*L9',L10,'s');
lmiterm([i,1,1,Q2],h2*L7',L12,'s');
lmiterm([i,1,1,N1],h2*(1/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],h2*(1/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],h2*(1/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],h2*(-1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],h2*(-1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],h2*(-1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,hatR2],h2*(1/h12)*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],h2*(-1/h12)*gamma3'*Coe,Coe*gamma3);

lmiterm([i,1,1,P],L1',L4,'s');
lmiterm([i,1,1,Q2],L8',L10,'s');
lmiterm([i,1,1,N1],(-h2/h12)*gamma2',gamma2,'s');
lmiterm([i,1,1,N2],(-h2/h12)*gamma2',gamma3,'s');
lmiterm([i,1,1,N3],(-h2/h12)*gamma2',gamma0,'s');
lmiterm([i,1,1,N4],(h1/h12)*gamma3',gamma2,'s');
lmiterm([i,1,1,N5],(h1/h12)*gamma3',gamma3,'s');
lmiterm([i,1,1,N3],(h1/h12)*gamma3'*(-1),gamma0,'s');
lmiterm([i,1,1,Q1],[e1',es'],[e1',es']');
lmiterm([i,1,1,Q1],-[e2',e11'],[e2',e11']');
lmiterm([i,1,1,Q2],L6',L6);
lmiterm([i,1,1,Q2],-L7',L7);
lmiterm([i,1,1,R1],es'*h1*h1,es);
lmiterm([i,1,1,R2],es'*h12*h12,es);
lmiterm([i,1,1,hatR1],-gamma1'*Coe,Coe*gamma1);
lmiterm([i,1,1,hatR2],-(2+(h1/h12))*gamma2'*Coe,Coe*gamma2);
lmiterm([i,1,1,hatR2],-(1-(h1/h12))*gamma3'*Coe,Coe*gamma3);


lmiterm([i,1,2,-N1],gamma2',1);
lmiterm([i,1,2,-N2],gamma3',1);
lmiterm([i,1,2,-N3],gamma0',1);

lmiterm([i,2,2,hatR2],-Coe,Coe);%t>=2



%%
i=i+1;

lmiterm([-i,1,1,P],1,1);

i=i+1;
lmiterm([-i,1,1,Q1],1,1);

i=i+1;
lmiterm([-i,1,1,Q2],1,1);

i=i+1;
lmiterm([-i,1,1,R1],1,1);

i=i+1;
lmiterm([-i,1,1,R2],1,1);
%% 求解语句

lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys);
if tmin < 0
    ACC1 = h2;
    if ACC2 - ACC1 <= acc
        disp (ACC1);
        break;
    end
else 
    ACC2 = h2;

 end

end
toc
