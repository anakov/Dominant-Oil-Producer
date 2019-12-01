set(0,'DefaultFigureWindowStyle','docked') 
global M_ 
% THE LINES BELOW ARE FOR THE NATURAL ALLOCATION - GENERATED AS THE ONE
% WHICH COINCIDES WITH STRICT INFLATION TARGETING
M_.params = repmat(NaN,42, 1);
M_.params(15) = 1;
phi_i = M_.params(15);
M_.params(12) = 5;   
phi_p_LR = M_.params(12);
M_.params(14) = 0;
phi_y_LR = M_.params(14);
M_.params(41) = 0;
phi_yg_LR = M_.params(41);
M_.params(13) = 0;
phi_o_LR = M_.params(13);
M_.params(42) = 0;
phi_oc_LR = M_.params(42);
M_.params(16) = 0;
phi_R = M_.params(16);

EOil4fig

nn = size(po_ea);
poem = 1;
pom  = oo_.mean(1);
Oem  = oo_.mean(2);
Om   = oo_.mean(3);
Xem  = oo_.mean(4);
Xm   = oo_.mean(5);
SO   = oo_.mean(6);
Yem  = oo_.mean(7);
Ym   = oo_.mean(8);

% IN THE EFFICIENT EQUILIBRIUM:
  % The price of oil responds only to Z:
poe_ez    = -zz_ez;      
poe_ea    = zeros(nn); % not to A 
poe_ex    = zeros(nn); % nor to X 
  % X doesn't affect Y nor O 
Ye_ex = zeros(nn); 
Oe_ex = zeros(nn); 
  % X only depends on X shock
Xe_ea = zeros(nn); 
Xe_ez = zeros(nn); 

% SOe is not computed inside
SOe_ea = (Oem + Oe_ea)./(Oem + Oe_ea + Xem + Xe_ea) - Oem./(Oem + Xem);
SOe_ez = (Oem + Oe_ez)./(Oem + Oe_ez + Xem + Xe_ez) - Oem./(Oem + Xem);
SOe_ex = (Oem + Oe_ex)./(Oem + Oe_ex + Xem + Xe_ex) - Oem./(Oem + Xem);


% PLOTS OF IMPULSE-RESPONSE FUNCTIONS

% AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
% figure(1)
% subplot(4,2,1);
% plot(100*a_ea,'b-o','MarkerSize',3);
% title('US technology shock (%)');
% axis tight
% 
% subplot(4,2,2);
% plot(100*poe_ea/poem,'k-^','Markersize',3) 
% hold on
% plot(100*po_ea/pom,'r-s','Markersize',3) 
% 
% subplot(4,2,3);
% plot(100*Oe_ea/Oem,'k-^','Markersize',3)
% hold on
% plot(100*O_ea/Om,'r-s','Markersize',3)
% 
% subplot(4,2,4);
% plot(100*Xe_ea/Xem,'k-^','Markersize',3)
% hold on
% plot(100*XX_ea/Xm,'r-s','MarkerSize',3)
% 
% subplot(4,2,5)
% plot(100*SOe_ea,'k-^','MarkerSize',3)
% hold on
% plot(100*SO_ea,'r-s','MarkerSize',3)
% 
% subplot(4,2,7);
% plot(100*Ye_ea/Yem,'k-^','Markersize',3)
% hold on
% plot(100*Y_ea/Ym,'r-s','Markersize',3)
% 
% subplot(4,2,8);
% %plot(100*(Ysgap_ea),'b-o','Markersize',3)
% hold on
% %plot(100*(Yngap_ea),'r-s','Markersize',3)


% ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

figure(2)
subplot(4,2,1);
plot(100*zz_ez,'b-o','MarkerSize',3);
title('Oil technology shock (%)');
axis tight

subplot(4,2,2);
plot(100*poe_ez/poem,'k-^','Markersize',3) 
hold on
plot(100*po_ez/pom,'r-s','Markersize',3) 

subplot(4,2,3);
plot(100*Oe_ez/Oem,'k-^','Markersize',3)
hold on
plot(100*O_ez/Om,'r-s','MarkerSize',3)

subplot(4,2,4);
plot(100*Xe_ez/Xem,'k-^','Markersize',3)
hold on
plot(100*XX_ez/Xm,'r-s','MarkerSize',3)

subplot(4,2,5)
plot(100*SOe_ez,'k-^','MarkerSize',3)
hold on
plot(100*SO_ez,'r-s','MarkerSize',3)

subplot(4,2,7);
plot(100*Ye_ez/Yem,'k-^','Markersize',3)
hold on
plot(100*Y_ez/Ym,'r-s','Markersize',3)

subplot(4,2,8);
plot(100*(Ygap_ez),'r-s','Markersize',3)
hold on

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

figure(3)
subplot(4,2,1);
plot(100*x_ex,'b-o','MarkerSize',3);
title('Non-OPEC/China shock (%)');
axis tight

subplot(4,2,2);
plot(100*poe_ex/poem,'k-^','Markersize',3) 
hold on
plot(100*po_ex/pom,'r-s','Markersize',3) 

subplot(4,2,3);
plot(100*Oe_ex/Oem,'k-^','Markersize',3)
hold on
plot(100*O_ex/Om,'r-s','MarkerSize',3)

subplot(4,2,4);
plot(100*Xe_ex/Xem,'k-^','Markersize',3)
hold on
plot(100*XX_ex/Xm,'r-s','MarkerSize',3)

subplot(4,2,5)
plot(100*SOe_ex,'k-^','MarkerSize',3)
hold on
plot(100*SO_ex,'r-s','MarkerSize',3)

subplot(4,2,7);
plot(100*Ye_ex/Yem,'k-^','Markersize',3)
hold on
plot(100*Y_ex/Ym,'r-s','Markersize',3)

subplot(4,2,8);
plot(100*(Ygap_ex),'r-s','Markersize',3)
hold on

figure(4)
subplot(5,2,3);
plot(400*PI_ez,'r-s','MarkerSize',3)
hold on

subplot(5,2,5);
plot(400*R_ez,'r-s','MarkerSize',3)
hold on

subplot(5,2,4);
plot(400*PI_ex,'r-s','MarkerSize',3)
hold on

subplot(5,2,6);
plot(400*R_ex,'r-s','MarkerSize',3)
hold on

subplot(5,2,1);
plot(100*zz_ez,'b-o','MarkerSize',3);
title('Oil technology shock (%)');
axis tight
subplot(5,2,7);
plot(100*Ye_ez/Yem,'k-^','Markersize',3)
hold on
plot(100*Y_ez/Ym,'r-s','Markersize',3)
subplot(5,2,9);
plot(100*(Ygap_ez),'r-s','Markersize',3)
hold on
subplot(5,2,2);
plot(100*x_ex,'b-o','MarkerSize',3);
title('Non-OPEC/China shock (%)');
axis tight
subplot(5,2,8);
plot(100*Ye_ex/Yem,'k-^','Markersize',3)
hold on
plot(100*Y_ex/Ym,'r-s','Markersize',3)
subplot(5,2,10);
plot(100*(Ygap_ex),'r-s','Markersize',3)
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THIS IS THE RULE OF INTEREST
M_.params( 15 ) = 1;
phi_i = M_.params( 15 );
M_.params( 12 ) = 2;
phi_p_LR = M_.params( 12 );
M_.params( 14 ) = 0;
phi_y_LR = M_.params( 14 );
M_.params( 13 ) = 0;
phi_o_LR = M_.params( 13 );
M_.params( 16 ) = 0.8;
phi_R = M_.params( 16 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EOil4fig

% AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
% figure(1)
% 
% subplot(4,2,2);
% plot(100*po_ea/pom,'b-o','MarkerSize',3)
% legend('po^e','po^n','po')
% title('Oil price (%)')
% 
% subplot(4,2,3);
% plot(100*O_ea/Om,'b-o','MarkerSize',3)
% legend('O^e','O^n','O')
% title('OPEC supply (%)')
% 
% subplot(4,2,4);
% plot(100*XX_ea/Xm,'b-o','MarkerSize',3)
% legend('X^e','X^n','X')
% title('Non-OPEC supply (%)')
% 
% subplot(4,2,5)
% plot(100*SO_ea,'b-o','MarkerSize',3)
% legend('s^e','s^n','s')
% title('OPEC share (pp)')
% 
% subplot(4,2,6)
% plot(400*PI_ea,'b-o','MarkerSize',3)
% hold on
% plot(400*R_ea,'k-s','MarkerSize',3)
% legend('PI','R',4)
% title('Inflation and interest rate (annualized pp)')
% 
% subplot(4,2,7);
% plot(100*Y_ea/Ym,'b-o','MarkerSize',3)
% legend('Y^e','Y^n','Y')
% title('US output (%)')
% 
% subplot(4,2,8);
% plot(100*(Ygap_ea),'k-^','MarkerSize',3)
% %legend('Y^{s}_{gap}','Y^{n}_{gap}','Y_{gap}',4)
% title('US output gap (pp)')


% ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ

figure(2)

subplot(4,2,2);
plot(100*po_ez/pom,'b-o','MarkerSize',3)
legend('po^e','po^n','po')
title('Oil price (%)')

subplot(4,2,3);
plot(100*O_ez/Om,'b-o','MarkerSize',3)
legend('O^e','O^n','O')
title('OPEC supply (%)')

subplot(4,2,4);
plot(100*XX_ez/Xm,'b-o','MarkerSize',3)
legend('X^e','X^n','X')
title('Non-OPEC supply (%)')

subplot(4,2,5)
plot(100*SO_ez,'b-o','MarkerSize',3)
legend('s^e','s^n','s')
title('OPEC share (pp)')

subplot(4,2,6)
plot(400*PI_ez,'b-o','MarkerSize',3)
hold on
plot(400*R_ez,'k-s','MarkerSize',3)
legend('PI','R',4)
title('Inflation and interest rate (annualized pp)')

subplot(4,2,7);
plot(100*Y_ez/Ym,'b-o','MarkerSize',3)
legend('Y^e','Y^n','Y')
title('US output (%)')

subplot(4,2,8);
plot(100*(Ygap_ez),'b-o','MarkerSize',3)
legend('Y^{n}_{gap}','Y_{gap}',4)
title('US output gap (pp)')

% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

figure(3)

subplot(4,2,2);
plot(100*po_ex/pom,'b-o','MarkerSize',3)
legend('po^e','po^n','po')
title('Oil price (%)')

subplot(4,2,3);
plot(100*O_ex/Om,'b-o','MarkerSize',3)
legend('O^e','O^n','O')
title('OPEC supply (%)')

subplot(4,2,4);
plot(100*XX_ex/Xm,'b-o','MarkerSize',3)
legend('X^e','X^n','X')
title('Non-OPEC supply (%)')

subplot(4,2,5)
plot(100*SO_ex,'b-o','MarkerSize',3)
legend('s^e','s^n','s')
title('OPEC share (pp)')

subplot(4,2,6)
plot(400*PI_ex,'b-o','MarkerSize',3)
hold on
plot(400*R_ex,'k-s','MarkerSize',3)
legend('PI','R',4)
title('Inflation and interest rate (annualized pp)')

subplot(4,2,7);
plot(100*Y_ex/Ym,'b-o','MarkerSize',3)
legend('Y^e','Y^n','Y')
title('US output (%)')

subplot(4,2,8);
plot(100*(Ygap_ex),'b-o','MarkerSize',3)
legend('Y^{n}_{gap}','Y_{gap}',4)
title('US output gap (pp)')


figure(4)

subplot(5,2,3);
plot(400*PI_ez,'b-o','MarkerSize',3)
legend('PI^{SIT}','PI^{FIT}',4)
title('Inflation (annualized pp)')

subplot(5,2,5);
plot(400*R_ez,'b-o','MarkerSize',3)
legend('R^{SIT}','R^{FIT}',4)
title('Interest rate (annualized pp)')

subplot(5,2,4);
plot(400*PI_ex,'b-o','MarkerSize',3)
legend('PI^{SIT}','PI^{FIT}',4)
title('Inflation (annualized pp)')

subplot(5,2,6);
plot(400*R_ex,'b-o','MarkerSize',3)
legend('R^{SIT}','R^{FIT}',4)
title('Interest rate (annualized pp)')

subplot(5,2,7);
plot(100*Y_ez/Ym,'b-o','MarkerSize',3)
legend('Y^e','Y^{SIT}','Y^{FIT}')
title('US output (%)')

subplot(5,2,9);
plot(100*(Ygap_ez),'b-o','MarkerSize',3)
legend('Y^{SIT}_{gap}','Y^{FIT}_{gap}',4)
title('US output gap (pp)')


subplot(5,2,8);
plot(100*Y_ex/Ym,'b-o','MarkerSize',3)
legend('Y^e','Y^{SIT}','Y^{FIT}')
title('US output (%)')

subplot(5,2,10);
plot(100*(Ygap_ex),'b-o','MarkerSize',3)
legend('Y^{SIT}_{gap}','Y^{FIT}_{gap}',4)
title('US output gap (pp)')

