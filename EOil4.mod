
var V_IMP, V_EXP, pomu, profit, Ri, pie, g_VA, g_po, SO, XX, Y, C, L, O, D, N, po, w, PI, R, DP, mc, 
    Ygap, Ye, Ce, Cgap, Oe, Xe, poe, L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, a, zz, x, rr, zb;  

varexo ea, ez, ex, er, ezb;

parameters A, Z, aa, bb, alpha, gamma, beta, theta, epsilon, mu, psi,
           phi_p_LR, phi_o_LR, phi_y_LR, phi_i, phi_R, R_bar, PI_bar, 
           rho_a, rho_z, rho_zb, rho_x, rho_r, stderrea, stderrezb, stderrez, stderrex, stderrer,
           po_ss, pomu_ss, Cgap_ss, Ygap_ss, SO_ss, X_O, K_ss, Lf_ss, Of_ss, Xf_ss, pof_ss, Omeg_fss, phi_yg_LR,
           phi_oc_LR;

% Rule switch: 0-Targeting; 1-Taylor 
% Note: interest rate targeting is when phi_i=1 and all other phi's are zero
phi_i =	1;   

% Taylor rule reaction coefficients (when phi_i=1) / target weights (when phi_i=0)
phi_p_LR  =	5;  % long-run inflation reaction / inflation weight
phi_y_LR  = 0;  % long-run output gap change reaction / output gap change weight
phi_yg_LR = 0;  % long-run output gap reaction / output gap weight
phi_o_LR  = 0.005;  % long-run oil price level reaction / oil price level weight
phi_oc_LR = 0;  % long-run oil price change reaction / oil price change weight

% Interest rate smoothing parameter
phi_R	  =	0;

beta      = 1.03^-.25;
PI_bar    = 1;
R_bar     = PI_bar/beta;

psi       =	1;
gamma     = 0.32; 
alpha     = 1-gamma-0.05;

theta	  =	0.75;
epsilon   = 7.66;
mu        = epsilon/(epsilon-1);

A         = 1; 
Z         = 1;
K_ss      = 1;
X_O       = 3/2;
bb        = 1;
aa        = 0;

Lf_ss = (alpha/(mu-(1-alpha-gamma)))^(1/(1+psi));  
pof_ss = (1+X_O*(1+alpha+gamma))/((1-alpha-gamma)+ X_O*(1+alpha+gamma));   
Of_ss  = (pof_ss*mu/((1-alpha-gamma)*Lf_ss))^(-1/(alpha+gamma))/(1+X_O);
Omeg_fss = X_O*Of_ss/(bb*(pof_ss-aa));

rho_a     = 0.95;
rho_z     = 0.70;
rho_x     = 0.80;
rho_zb	  =	0;
rho_r     = 0;

stderrea  = 0.01;
stderrez  = 0.15;
stderrex  = 0.30;
stderrezb =	0;
stderrer  = 0;        % 0.0025/4*1.65658;  %  shut off monetary shock for welfare evaluation


% analytic expression for flex-price steady-state labor
Lf_ss = (alpha/(mu-(1-alpha-gamma)))^(1/(1+psi));  

% pre-solving the steady-state
X_ss = fsolve('EOil_do_SS', [1.45, 0],[],'zerofinder');
[X_ss,SScheck]=EOil_do_SS(X_ss,'monopolist');              

% retrieving necessary steady-state values
pomu_ss = X_ss(3);
SO_ss   = X_ss(9);
po_ss   = X_ss(17);
Cgap_ss = X_ss(26);

model;
# phi_p     = phi_p_LR*(1-phi_R);    % short-run inflation coefficient
# phi_y     = phi_y_LR*(1-phi_R);    % short-run output gap change coefficient
# phi_yg    = phi_yg_LR*(1-phi_R);   % short-run output gap coefficient
# phi_o     = phi_o_LR*(1-phi_R);    % short-run oil price level coefficient
# phi_oc    = phi_oc_LR*(1-phi_R);   % short-run oil price change coefficient

V_IMP = log(C) - L^(1+psi)/(1+psi) + beta*V_IMP(+1);

V_EXP  = log(O*(po-poe)) + beta*V_EXP(+1);

pomu = po/poe;

SO = O/(XX+O);

XX = bb*(po*exp(zz)-aa)*exp(x)*Omeg_fss;

C = Y - po*(O + XX);

w = C*L^psi;

1 = beta*R*exp(zb(+1))*C/C(+1)/exp(zb)/PI(+1);

D = exp(zb)*Y/C + beta*theta*PI_bar*PI(+1)^(epsilon-1)*D(+1);

N = mu*mc*exp(zb)*Y/C + beta*theta*PI(+1)^epsilon*N(+1);

1 = theta*(PI/PI_bar)^(epsilon-1) + (1-theta)*(N/D)^(1-epsilon);

DP = theta*(PI/PI_bar)^epsilon*DP(-1) + (1-theta)*(N/D)^(-epsilon);

po*(O + XX) = (1-alpha-gamma)*mc*Y*DP;

w*L = alpha*mc*Y*DP;

Y*DP = exp(a)*L^alpha*K_ss^gamma*(O + XX)^(1-alpha-gamma);
                                                                                                             
phi_i*(R-R_bar)= phi_i*(R(-1)-R_bar)*phi_R + (PI-1)*phi_p + (C-Ce-(C(-1)-Ce(-1)))*phi_y + 
                 + (C/Ce-Cgap_ss)*phi_yg + (po-po_ss)*phi_o + (po-po(-1))*phi_oc +rr; 

0 = 1/O  - L1*po - L7*po + L9*(1-alpha-gamma)*Y*DP/(O + XX);

0 = - L1*C + L2 - L2(-1)*R(-1)*C(-1)/C/PI - L3*Y -L4*mu*mc*Y + L10*w + (L11-beta*L11(+1))*phi_y*C
      + L11*phi_yg*C/Ce;

0 =  L1 + L3*exp(zb) + L4*mu*mc*exp(zb) + L7*(1-alpha-gamma)*mc*DP - L8*alpha*mc*DP - L9*DP;

0 = L3(-1)*theta*PI_bar*C(-1)*PI^(epsilon-1) - L3*C + L5*(1-theta)*(epsilon-1)*N^(1-epsilon)*D^(epsilon-2) 
    + L6*(1-theta)*epsilon*N^(-epsilon)*D^(epsilon-1);

0 = L4(-1)*theta*C(-1)*PI^(epsilon) - L4*C + L5*(1-theta)*(1-epsilon)*N^(-epsilon)*D^(epsilon-1)
    - L6*(1-theta)*epsilon*N^(-epsilon-1)*D^epsilon;

0 = L8*w + L9*alpha*Y*DP/L + L10*C*psi*L^(psi-1);

0 = 1/(po - 1/exp(zz)) - (O + 2*XX)*(L1 + L7) + L9*(1-alpha-gamma)*Y*DP/(O + XX)*exp(zz)*exp(x)*Omeg_fss*bb
    +  phi_o*L11 + phi_oc*(L11-beta*L11(+1)) ; 

0 = -L2(-1)*R(-1)*C(-1)/C/PI^2 + L3(-1)*theta*PI_bar*(epsilon-1)*C(-1)*PI^(epsilon-2)*D 
   + L4(-1)*theta*epsilon*C(-1)*PI^(epsilon-1)*N + L5*theta*(epsilon-1)*PI^(epsilon-2)*PI_bar^(1-epsilon) 
   + L6*theta*epsilon*PI_bar^(-epsilon)*PI^(epsilon-1)*DP(-1)
   + L11*phi_p;

0 = L6(+1)*beta*theta*(PI(+1)/PI_bar)^epsilon - L6 + L7*(1-alpha-gamma)*mc*Y - L8*alpha*mc*Y - L9*Y;

0 = L4*mu*exp(zb) + L7*(1-alpha-gamma)*DP - L8*alpha*DP;

0 = L2/R - phi_i*(L11 - phi_R*beta*L11(+1));

0 = L8*L - L10;

Xe = bb*(poe*exp(zz)-aa)*exp(x)*Omeg_fss;

poe*(Oe + Xe) = (1-alpha-gamma)/mu*Ye;

poe = 1/exp(zz);

Ye = exp(a)*Lf_ss^alpha*K_ss^gamma*(Oe+Xe)^(1-alpha-gamma);

Ce = Ye*(1-(1-alpha-gamma)/mu);

Ygap = Y/Ye;

Cgap = C/Ce;

a = rho_a*a(-1) + ea;

zz = rho_z*zz(-1) - ez;

x = rho_x*x(-1) - ex;

rr = rho_r*rr(-1) - er;

zb = rho_zb*zb(-1) + ezb;

g_VA = C/C(-1) - 1;

pie = PI - 1;

Ri = R - 1/beta;

g_po = po/po(-1) - 1;

profit = O*(po-poe);

end;
initval;
    V_IMP     =	    X_ss(1);
    V_EXP     =	    X_ss(2);
    pomu      =     X_ss(3);
    profit      = 	X_ss(4);
    Ri        = 	X_ss(5);
    pie       = 	X_ss(6);
    g_VA      = 	X_ss(7);
    g_po      = 	X_ss(8);
    SO        = 	X_ss(9);
    XX        = 	X_ss(10);
    Y         = 	X_ss(11);
    C         = 	X_ss(12);
    L         = 	X_ss(13);
    O         = 	X_ss(14);
    D         = 	X_ss(15);
    N         = 	X_ss(16);
    po        = 	X_ss(17);
    w         = 	X_ss(18);
    PI        = 	X_ss(19);
    R         = 	X_ss(20);
    DP        = 	X_ss(21);
    mc        = 	X_ss(22);
    Ygap      = 	X_ss(23);
    Ye        = 	X_ss(24);
    Ce        = 	X_ss(25);
    Cgap      = 	X_ss(26);
    Oe        = 	X_ss(27);
    Xe        = 	X_ss(28);
    poe       = 	X_ss(29);
    L1        = 	X_ss(30);
    L2        = 	X_ss(31);
    L3        = 	X_ss(32);
    L4        = 	X_ss(33);
    L5        = 	X_ss(34);
    L6        = 	X_ss(35);
    L7        = 	X_ss(36);
    L8        = 	X_ss(37);
    L9        = 	X_ss(38);
    L10       = 	X_ss(39);
    L11       = 	X_ss(40);
    a         = 	X_ss(41);
    zz        = 	X_ss(42);
    x         = 	X_ss(43);
    rr        = 	X_ss(44);
    zb        = 	X_ss(45);
end;

steady;
% check

shocks;
var ea  = stderrea^2;
var ez  = stderrez^2;
var ezb = stderrezb^2;
var ex  = stderrex^2;
var er  = stderrer^2;
end;

if 0
stoch_simul(order=2, nograph,irf=0) g_VA pie Ri g_po SO Cgap po O profit XX V_IMP V_EXP;
stddevs=sqrt(diag(oo_.var));
stddevsvec=[stddevs(1:5)];
corr1=oo_.var(1,4)/sqrt(oo_.var(1,1)*oo_.var(4,4));
corr2=oo_.var(2,3)/sqrt(oo_.var(2,2)*oo_.var(3,3));
corr3=oo_.var(5,4)/sqrt(oo_.var(5,5)*oo_.var(4,4));
corr4=oo_.var(3,4)/sqrt(oo_.var(3,3)*oo_.var(4,4));
corr5=oo_.var(3,1)/sqrt(oo_.var(3,3)*oo_.var(1,1));
corr6=oo_.var(2,1)/sqrt(oo_.var(2,2)*oo_.var(1,1));
corr7=oo_.var(4,9)/sqrt(oo_.var(4,4)*oo_.var(9,9));
corr8=oo_.var(8,10)/sqrt(oo_.var(8,8)*oo_.var(10,10));
corrvec=[corr1;corr2;corr3;corr4;corr5;corr6;corr7;corr8];
autocorrvec=diag(oo_.autocorr{1});
autocorrvec=autocorrvec(1:5);
momentsvec=[stddevsvec;NaN;corrvec;NaN;autocorrvec];
disp(sprintf('\n'));
disp(momentsvec);
end

% use next line for IRFs
  stoch_simul(order=1, nograph, irf=20) po Oe O Xe XX SO Ye Y Cgap rr a zz x Ygap pomu PI R g_VA;

% use next line for 2-nd order stochastic simulations
% stoch_simul(order=2, nograph, irf=0) pie Cgap C L profit V_IMP V_EXP g_po;

stddevs=(diag(oo_.var)).^.5;
% means=[oo_.mean(6);oo_.mean(7);400*oo_.mean(1); 100*oo_.mean(2)/Cgap_ss; oo_.mean(3); oo_.mean(4) ; 100*oo_.mean(5)/X_ss(4) ; oo_.mean(8) ; oo_.mean(9)]
% disp('Std devs')
% 100*[4*stddevs(1);stddevs(2);stddevs(3);stddevs(4);stddevs(5);stddevs(8);stddevs(9)]  

