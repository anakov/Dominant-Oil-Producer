% Computes steady-state (non-Dynare format):
% case 'zerofinder' 
  % input: [po, L11] - initial guess for price of oil and Lagrange multiplier 
  % output: steady-state values of [po, L11] 
% case 'monopolist'
  % input: steady-state values of [po, L11] 
  % output: rest of steady-state values 

function [stst_vec, SScheck] = EOil_do_SS(xini,market_case)
global M_

SScheck  = -inf;

% needed parameters only
gamma    = M_.params(6);
beta     = M_.params(7);
epsilon  = M_.params(9);
mu       = M_.params(10);
PI_bar   = M_.params(18);
Omeg_fss = M_.params(40);
A        = M_.params(1);
Z        = M_.params(2);
K_ss     = M_.params(35);
bb       = M_.params(4);
aa       = M_.params(3);
alpha    = M_.params(5);
theta    = M_.params(8);
psi      = M_.params(11);

% Steady state values independent of oil price
PI  = PI_bar;
N_D = ((1-theta*PI^(epsilon-1))/(1-theta))^(1/(1-epsilon));
mc  = 1/mu*N_D*(1-beta*theta*PI^(epsilon))/(1-beta*theta*PI^(epsilon-1));
DP  = (1-theta)*N_D^(-epsilon)/(1-theta*PI^(epsilon));
L   = (alpha*mc*DP/(1-(1-alpha-gamma)*mc*DP))^(1/(1+psi));
R   = PI/beta;

switch market_case
    case{'zerofinder'}  % used to pre-solve steady-state
        po  = xini(1);
        % Setting L11=0 we are imposing that the monopolist ignores the central bank's threat.
        % This also makes it easier to solve for this (restricted) steady-state
        % when the MH algorithm varies the parameters psi, theta and phi_R.
        % Could speed up slightly the steady-state computation by replacing L11=0 in
        % equations for L2, L3, etc
        L11 = 0 ;
           % L11 = xini(2);
        Od = ((1-alpha-gamma)*mc*A*L^alpha*K_ss^gamma/po)^(1/(alpha+gamma));
        Y  = A*L^alpha*K_ss^gamma*Od^(1-alpha-gamma)/DP;
        C  = Y*(1-(1-alpha-gamma)*mc*DP);
        XX = Omeg_fss*po*Z;
        O  = Od-XX;
        D  = Y/C/(1-beta*theta*PI^(epsilon-1));
        N  = N_D*D;
       
        if (O>=0 & po>1/Z+1e-9)
            L9  = mc*(1/(po*Z-1)-2*XX/O)/po/Od;
            L8  = -L9/mc/(1+psi);
            L10 = L8*L;
            L2  = 0;
            L3  = (-L2*(1-1/beta)/C-L10*L^psi+L9*DP)/(1-Y/C*(1-DP*mu*mc));
            L4  = -L3*DP;
            L1  = -L3*Y/C*(1-mu*mc*DP)+L2*(1-1/beta)/C+L10*L^psi;
            L7  = (L8*alpha+L3*mu)/(1-alpha-gamma);
            L6  = ((1-alpha-gamma)*L7*mc*Y-L8*alpha*mc*Y-L9*Y)/(1-beta*theta*PI^epsilon);
            L5  = D/(epsilon-1)*(L3*C-L6*epsilon/N);
        
          % function output
            stst_vec(1) = 1/O-(L1+L7)*po+L9*po/mc;
            stst_vec(2) = L3*C*theta*(epsilon-1)*PI^(epsilon-1)*D+L4*theta*epsilon*C*N*PI^epsilon+...
                          L5*(epsilon-1)*theta*PI^(epsilon-1)+L6*theta*epsilon*DP*PI^epsilon;
            SScheck = 0;
            
        elseif (O<0 & po>1/Z+1e-9)
          % function output
            stst_vec = po-1/Z;
        else
          % function output 
            stst_vec = 1e16*(po-1/Z+1e-9)*(po>=1/Z); %penalization
        end
        
    case{'monopolist'}
        po  = xini(1);
        % setting L11=0 we are imposing that the monopolist ignores the central bank's threat
        % this also makes it easier to solve for this (restricted) steady-state
        % when the MH algorithm varies the parameters psi, theta and phi_R
        % Could speed up slightly the steady-state computation by replacing L11=0 in
        % equations for L2, L3, etc
        L11 = 0; 
            % L11 = xini(2); 
        Od  = ((1-alpha-gamma)*mc*A*L^alpha*K_ss^gamma/po)^(1/(alpha+gamma));
        Y   = A*L^alpha*K_ss^gamma*Od^(1-alpha-gamma)/DP;
        C   = Y*(1-(1-alpha-gamma)*mc*DP);
        w   = C*L^psi;
        XX  = Omeg_fss*po*Z;
        O   = Od-XX;
        D   = Y/C/(1-beta*theta*PI^(epsilon-1));
        N   = N_D*D;
        
        if (O>0 & po>1/Z+1e-9)
            L9  = mc*(1/(po*Z-1) - 2*XX/O)/po/Od ; 
            L8  = -L9/mc/(1+psi);
            L10 = L8*L;
            L2  = 0;
            L3  = (-L2*(1-1/beta)/C-L10*L^psi+L9*DP)/(1-Y/C*(1-DP*mu*mc));
            L4  = -L3*DP;
            L1  = L9*DP - L3;
            L7  = (L8*alpha+L3*mu)/(1-alpha-gamma);
            L6  = ((1-alpha-gamma)*L7*mc*Y-L8*alpha*mc*Y-L9*Y)/(1-beta*theta*PI^epsilon);
            L5  = D/(epsilon-1)*(L3*C-L6*epsilon/N);
            SScheck = 0;
        else
            po = ((1-alpha-gamma)*mc*A*L^alpha*K_ss^gamma/(Omeg_fss*Z)^(alpha+gamma))^(1/(1+alpha+gamma));
            XX = Omeg_fss*po*Z;
            Od = XX;
            Y  = A*L^alpha*K_ss^gamma*Od^(1-alpha-gamma)/DP;
            C  = Y*(1-(1-alpha-gamma)*mc*DP);
            w  = C*L^psi;
            O  = 0;
            L1 = 0;
            L2 = 0;
            L3 = 0;
            L4 = 0;
            L5 = 0;
            L6 = 0;
            L7 = 0;
            L8 = 0;
            L9 = 0;
            L10 = 0;
            L11 = 0;
        end
        SO  = O/Od;
        a  = 0;
        rr = 0;
        x  = 0;
        zz = 0;
        zb = 0;
        g_VA = 0;
        pie  = 0;
        Ri   = 0;
        g_po = 0;
      
        Ye   = (L^alpha*K_ss^gamma*((1-alpha-gamma)/mu)^(1-alpha-gamma))^(1/(alpha+gamma));
        Ce   = Ye*(1-(1-alpha-gamma)/mu);
        Ygap = Y/Ye;
        Cgap = C/Ce;
        poe  = 1/Z;
        Xe   = poe*Omeg_fss;
        Oe   = ((1-alpha-gamma)/mu*L^alpha*K_ss^gamma)^(1/(alpha+gamma))-Xe;
        pomu = po/poe;
       
        profit = O*(po-poe);
        V_IMP = (log(C) - L^(1+psi)/(1+psi))/(1-beta);
        V_EXP  = log(profit)/(1-beta);

      % function output
        stst_vec = [V_IMP, V_EXP, pomu, profit, Ri, pie, g_VA, g_po, SO, XX,...
                   Y, C, L, O, D, N, po, w, PI, R,...
                   DP, mc, Ygap, Ye, Ce, Cgap, Oe, Xe, poe, L1,...
                   L2, L3, L4, L5, L6, L7, L8, L9, L10, L11,...
                   a, zz, x, rr, zb]';
end

