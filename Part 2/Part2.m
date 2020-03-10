P_tan = date_indiv(218)
omeg = logspace(-2,2,1000)';

%  nyquist(P_tan,omeg);
[re,im] = nyquist(P_tan,omeg);
re = squeeze(re);
im = squeeze(im);

inters1 = P_tan.num{1}(4) / P_tan.den{1}(4) ;

for i = 2 : length(omeg)
    if  sign(im(i)) ~= sign(im(i-1)) 
        if abs(im(i)) < abs(im(i-1)) 
            inters2 = re(i);
        else
            inters2 = re(i-1);
        end
    end
end

%  nyquist(P_tan * 2 , omeg);
[re,im] = nyquist(P_tan * 2 , omeg) ;

re = squeeze(re);
im = squeeze(im);
for i = 2 : length(omeg)
    if sign(im(i)) ~= sign(im(i-1))
        if abs(im(i)) < abs(im(i-1))
            inters3 = re(i) ;
        else
            inters3 = re(i-1) ;
        end
    end
end
hold
% nyquist(P_tan,omeg);

%  nyquist(P_tan * exp(-1i * pi / 4 ) , omeg );
[re ,im ] = nyquist(P_tan * exp( -1i * pi / 4 ), omeg );

re = squeeze(re);
im = squeeze(im);

for i = 2 : length(omeg)
    if sign(im(i)) ~= sign(im(i-1))
        if abs(im(i)) < abs(im(i-1))
            inters4 = re(i) ;
        else
            inters4 = re(i-1) ;
        end
    end
end

hold on 
% nyquist(P_tan , omeg); 

s = tf('s');

 %nyquist(P_tan * (1/s) , omeg );
% [re,im] = nyquist(P_tan*(1/s) , omeg);

asimpt = -0.0905;
            
%Ex 2

K_1 = ( 100 * P_tan.den{1}(4) ) / P_tan.num{1}(4) ;
T_1 = 250 ;
C_1 = tf( K_1, [T_1 1 ])

%   nyquist(P_tan * C_1,omeg)
[re,im] = nyquist(P_tan*C_1,omeg);
re = squeeze(re);
im = squeeze(im);
for i = 2 : length(omeg)
    if sign(im(i)) ~= sign(im(i-1))
        if abs(im(i)) < abs(im(i-1))
            verificare = re(i) ;
        else
            verificare = re(i-1) ;
        end
    end
end
K_2 = ( 100 * P_tan.den{1}(4) ) / P_tan.num{1}(4) ; ;
T_2 = 800 ;
C_2 = K_2 * tf([1,1] , [T_2,1])

% nyquist(P_tan * C_2 , omeg);
% hold on; 
% nyquist(tf([-1 1], [1 1]), omeg);
[re,im] = nyquist(P_tan*C_2,omeg);
re = squeeze(re);
im = squeeze(im);
for i = 2 : length(omeg)
    if sign(im(i)) ~= sign(im(i-1))
        if abs(im(i)) < abs(im(i-1))
            verificare1 = re(i) ;
        else
            verificare1 = re(i-1) ;
        end
    end
end

%Ex 3 
%a)


% URSS ISI AMINTESTE DE VARIANTELE ANTERIOARE ALE ACESTUI EXERCITIU.


% % bode(P_tan,omeg);
% [amp,def,~] = bode(P_tan,omeg);
% U = ((sqrt(2)*s)/2 + sqrt(2)/2)/(s^2 + 1 ) 
% %  bode(P_tan*U,omeg)
% amp1 = amp(500) * 7 ;
% def1 = def(500) + pi/4 ;
% amp1 = 10^(amp1/20);
% 
% bode(3*P_tan,omeg);
% [amp,def,~] = bode(3*P_tan,omeg);
% U = ((sqrt(2)*s)/2 + sqrt(2)/2)/(s^2 + 1 ) 
% %  bode(P_tan*U,omeg)
% amp2 = amp(500) * 7 ;
% def2 = def(500) + pi/4 ;
% amp2 = 10^(amp2/20);

%Ex 3 sub a)
 U = ((sqrt(2)*s)/2 + sqrt(2)/2)/(s^2 + 1 ) 
 [amp,def,~] = bode(P_tan*U,omeg);
%  bode(P_tan*U,omeg);
 amp1 = 7* abs(evalfr(P_tan,1i));
 def1 = def(500)
 

 %Ex 3 sub b)
 
  U = ((sqrt(2)*s)/2 + sqrt(2)/2)/(s^2 + 1 ) 
 [amp,def,~] = bode(3*P_tan*U,omeg);
%  bode(3*P_tan*U,omeg);
amp2 = 7 * abs(evalfr(P_tan * 3,1i));
 def2 = def(500)


%Ex 3 sub c)

  U = ((sqrt(2)*s)/2 + sqrt(2)/2)/(s^2 + 1 ) 
 [amp,def,~] = bode(exp(-1i*pi/6)*P_tan*U,omeg);
%  bode(P_tan*exp(-1i*pi/6)*U,omeg);
 amp3 = 7 * abs(evalfr(P_tan * exp(-1i*pi/6),1i));
 def3 = def(500)
 
%Ex 3 sub d)

% [amp,def,~] = bode(100*P_tan,omeg);
%  bode(100*P_tan,omeg)
omeg_1 = 3.49;

%Ex 3 sub e)

% [amp,def,~] = bode(100*P_tan,omeg);
omeg_2 = 1.99;

%Ex 4 sub a)
K_3 = K_1;
w_3 = 1/5;
C_3 = tf(K_3 * w_3, [1 w_3]);
%   bode(P_tan * C_3, omeg);
 
%Ex 4 sub b)
A_4 = 1;
B_4 = omeg_2^2;
C_4 = 100 * tf(B_4 * [1 A_4], A_4 * [1 B_4]);
%  bode(P_tan * C_4, omeg);    %Faza este -143

save('Chesnoiu_Alex-Marian_322AC_tema2.mat', 'inters1', 'inters2', ...
'inters3', 'inters4', 'asimpt', 'K_1', 'T_1', 'K_2', 'T_2', ...
'amp1', 'def1', 'amp2', 'def2', 'amp3', 'def3', 'omeg_1', ...
'omeg_2', 'K_3', 'w_3', 'A_4', 'B_4');


