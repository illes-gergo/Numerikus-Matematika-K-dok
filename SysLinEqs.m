clc; close all; clear;

A = [1,-1,3;2,-1,2;2,-1,1];
b0 = [0;-4;-7];

% LDU felbontás
[L,U,P] = lu(A);

D=diag(diag(U));

for ii = 1:size(U,2)
    D2(ii,ii) = U(ii,ii);
end

Unew=U./diag(U);
DU=D*Unew;

disp(table(U,D,D2,Unew,DU));

%% Cholesky beépített függvénnyel
clear;clc;
A=[4 12 -16; 12 37 -43;-16 -43 98];
L=chol(A,'lower');
disp(table(A,L,L',L*L','VariableNames',{'A','L','Lt','L*Lt'}))

%% Cholesky általánosan

clear;clc;
A=sym('A',[3,3])
A=tril(A)+tril(A,-1).'
L=tril(sym('L',[3,3]));
disp(table(A,L));

eq=A==L*transpose(L)

%% Cholesky algoritmussal
clear;
A=[4 12 -16; 12 37 -43; -16 -43 98];
L=zeros(size(A));
for k1=1:size(A,1) % loop for j
for k2=(k1+1):size(A,2) % loop for i
L(k1,k1)=sqrt(A(k1,k1)-sum(L(k1,1:(k1-1)).*conj(L(k1,1:(k1-1)))));
L(k2,k1)=1/L(k1,k1).*(A(k2,k1)-sum(L(k2,1:(k1-1)).*conj(L(k1,1:(k1-1)))));
end
end
disp(table(A,L,L',L*L','VariableNames',{'A','L','Lt','L*Lt'}))

%% LDL Általánosan
clear;clc;
A=sym('A',[3,3]);
A=tril(A)+tril(A,-1).'; % A szimmetrikus mátrix
L=tril(sym('L',[3,3]),-1)+sym(diag(ones(3,1)));
D=diag(diag(sym('D',[3,3])));
disp(table(A,L,D));

eq=A==L*D*transpose(L)

%% LDL algoritmus
clear;
A=[4 12 -16; 12 37 -43; -16 -43 98];
L=zeros(size(A))+diag(ones(3,1));
D=zeros(size(A));
for k1=1:size(A,1) % loop for j
for k2=(k1+1):size(A,2) % loop for i
D(k1,k1)=A(k1,k1)-sum(L(k1,1:(k1-1)).*L(k1,1:(k1-1)).*diag(D(1:(k1-1),1:(k1-1))));
L(k2,k1)=(1./D(k1,k1)).*(A(k2,k1)-sum(L(k2,1:(k1-1)).*L(k1,1:(k1-1)).*diag(D(1:(k1-1),1:(k1-1)))));
end
end
disp(table(A,L,D,L',L*D*L','VariableNames',{'A','L','D','Lt','L*D*Lt'}))
