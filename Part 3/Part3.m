clear;clc;
[~, P_gir] = date_indiv(218);
P_gir = tf(ss(P_gir,'min'));
%Exercitiul 1
a1 = 3 ;
[X,Y,N,M] = eucl_Youla(P_gir.num{1},P_gir.den{1},a1);
s = tf('s');
Q1 = 1 / (s + 1) ;
C = (X + M*Q1)/(Y-N*Q1);
C1 = tf(ss(C,'min'));
T1 = tf(ss(P_gir*C1/(1+P_gir*C1),'min'));
%step(T1)
%V1 = stepinfo(T1)

%Exercitiul 2

a2 = 3 ;
[X2,Y2,N2,M2] = eucl_Youla(P_gir.num{1},P_gir.den{1},a2);
Q0 = evalfr(Y2,0)/evalfr(N2,0);
Q2 = 1 / (s + Q0^(-1));
C2 = (X2 + M2*Q2) / ( Y2 - N2*Q2) ;
C2 = tf(ss(C2,'min'));
T2 = P_gir*C2 /(1+P_gir*C2);
T2 = tf(ss(T2,'min'));
%step(T2)
%V2 = stepinfo(T2)

%Exercitiul 3

a3 = 3;
[X3,Y3,N3,M3] = eucl_Youla(P_gir.num{1},P_gir.den{1},a3);
Q3 = 1 / ( s + 1 ) ;
C3 = (X3+M3*Q3) / ( Y3 - Q3*N3);
C3 = tf(ss(C3,'min'));
%margin(P_gir*C3)
T3 = P_gir * C3 / ( 1 + P_gir*C3);
T3 = tf(ss(T3,'min'));
%V3 = stepinfo(T3)
%Exercitiul 4

a4 = 3 ;
[X4,Y4,N4,M4] = eucl_Youla(P_gir.num{1},P_gir.den{1},a4);
Q0 = evalfr(Y4,0)/evalfr(N4,0);
Q4 = 1 / (s+Q0^(-1));
C4 = ( X4 + M4*Q4) / ( Y4 - N4*Q4);
C4 = tf(ss(C4,'min'));
margin(P_gir*C4);
T4 = P_gir * C4 / ( P_gir*C4 + 1 ) ;
V4 = stepinfo(T4)

save('Chesnoiu_Alex-Marian_322AC_tema3.mat', 'a1', 'Q1', ...
    'a2', 'Q2', 'a3', 'Q3', 'a4', 'Q4');

