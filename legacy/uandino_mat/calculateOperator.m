function y = calculateOperator(neutrinoEnergy, A, L)
kuoPantaleone = 0;
ohlsson = 0;
real=1;
if kuoPantaleone
dM32 =1e-4;
dm21 = 1e-8;
theta1 = asin(sqrt(0.437));
theta2 = asin(sqrt(0.0214));
theta3 = asin(sqrt(0.297));
%{
theta1 = deg2rad(45);
theta2 = deg2rad(5.);
theta3 = deg2rad(45.);
%}
elseif ohlsson
dM32 = 3.2e-3;
dm21 = 0.;
theta1 = deg2rad(45.);
theta2 = deg2rad(5.);
theta3 = deg2rad(45.);
elseif real
dM32 = 2.53685e-3;
dm21=7.37e-5;
theta1 = asin(sqrt(0.437));
theta2 = asin(sqrt(0.0214));
theta3 = asin(sqrt(0.297));
end
CP=0;
if CP
    CP_phase = pi*1.35;
else
    CP_phase =0;
end

%U matrix Elements
Ue1 = cos(theta2)*cos(theta3);
Ue2 = sin(theta3)*cos(theta2);
Ue3 = sin(theta2)*exp(-1j*CP_phase);
Umu1=-sin(theta3)*cos(theta1)-sin(theta1)*sin(theta2)*cos(theta3)*exp(-1j*CP_phase);
Umu2=cos(theta1)*cos(theta3)-sin(theta1)*sin(theta2)*sin(theta3)*exp(-1j*CP_phase);
Umu3=sin(theta1)*cos(theta2);
Ut1=sin(theta1)*sin(theta3)-sin(theta2)*cos(theta1)*cos(theta3)*exp(-1j*CP_phase);
Ut2=-sin(theta1)*cos(theta3)-sin(theta2)*sin(theta3)*cos(theta1)*exp(-1j*CP_phase);
Ut3=cos(theta1)*cos(theta2);
CKM = [Ue1, Ue2, Ue3;Umu1, Umu2, Umu3;Ut1, Ut2, Ut3];



E21 = dm21/(2*neutrinoEnergy);
E32 = dM32/(2*neutrinoEnergy);
E12=-E21;
E23=-E32;
E31=E12-E23;
E13=-E31;
%Elements of the Tmatrix in mass basis
T_11=A*Ue1*Ue1-(1./3)*A+(1./3)*(E12+E13);
T_12=A*Ue1*Ue2;
T_13=A*Ue1*Ue3;
T_21=T_12;
T_22=A*Ue2*Ue2-(1./3)*A+(1./3)*(E21+E23);
T_23=A*Ue2*Ue3;
T_31=T_13;
T_32=T_23;
T_33=A*Ue3*Ue3-(1./3)*A+(1./3)*(E31+E32);
T_mass_mat = [T_11, T_12, T_13;T_21, T_22, T_23;T_31, T_32, T_33];
T_sq_11=(1./3)*(A*A*(Ue1*Ue1+(1./3))+2*A*(Ue1*Ue1-(1./3))*(E12+E13)+(1./3)*(E12+E13)*(E12+E13));
T_sq_12=(1./3)*Ue1*Ue2*A*(A+E13+E23);
T_sq_13=(1./3)*Ue1*Ue3*A*(A+E12+E32);
T_sq_21=T_sq_12;
T_sq_22=(1./3)*(A*A*(Ue2*Ue2+(1./3))+2*A*(Ue2*Ue2-(1./3))*(E21+E23)+(1./3)*(E21+E23)*(E21+E23));
T_sq_23=(1./3)*Ue2*Ue3*A*(A+E21+E31);
T_sq_31=T_sq_13;
T_sq_32=T_sq_23;
T_sq_33=(1./3)*(A*A*(Ue3*Ue3+(1./3))+2*A*(Ue3*Ue3-(1./3))*(E31+E32)+(1./3)*(E31+E32)*(E31+E32));
T_sq_mass_mat = [T_sq_11, T_sq_12, T_sq_13;T_sq_21, T_sq_22, T_sq_23;T_sq_31, T_sq_32, T_sq_33];
T_flav_mat = CKM*T_mass_mat*CKM';
T_sq_flav_mat = CKM*T_sq_mass_mat*CKM';
%Calculate c's
c1 = -A*A/3 + (A/(6*neutrinoEnergy))*(Ue1*Ue1*(dM32+2*dm21)+Ue2*Ue2*(dM32-dm21)-Ue3*Ue3*(2*dM32+dm21))-(1./(12*neutrinoEnergy*neutrinoEnergy))*(dM32*dM32+dm21*dm21+dM32*dm21);
c0 = (-2./27)*A*A*A+(A*A/(18*neutrinoEnergy))*(Ue1*Ue1*(dM32+2*dm21) + Ue2*Ue2*(dM32-dm21) - Ue3*Ue3*(2*dM32+dm21))+(A/(36*neutrinoEnergy*neutrinoEnergy))*(Ue1*Ue1*(2*dM32+dm21)*(dM32-dm21)+Ue2*Ue2*(2*dM32+dm21)*(dM32+2*dm21)-Ue3*Ue3*(dM32+2*dm21)*(dM32-dm21)-(dM32*dM32+dm21*dm21+dM32*dm21))-(1./(216*neutrinoEnergy*neutrinoEnergy*neutrinoEnergy))*(2*dM32+dm21)*(dM32+2*dm21)*(dM32-dm21);

%Eigenvals

s1PlusS2 = 2*sqrt((-1./3)*c1)*cos((1./3)*atan((1./c0)*sqrt(-c0*c0-(4./27)*c1*c1*c1)));
s1MinusS2 = -2j*sqrt((-1./3)*c1)*sin((1./3)*atan((1./c0)*sqrt(-c0*c0-(4./27)*c1*c1*c1)));
lam1 = -0.5*s1PlusS2 + sqrt(3.)*1j*s1MinusS2/2;
lam2 = -0.5*s1PlusS2 - sqrt(3.)*1j*s1MinusS2/2;
lam3 = s1PlusS2;
lam = [lam1, lam2, lam3];

%Calculate operator

trace_hamiltonian=0.5*E21+E32+3*neutrinoEnergy+A;
phi_phase = exp(-1j*L*trace_hamiltonian);
summ = zeros(3,3);
for a=1:3
   summ = summ + exp(-1j*L*lam(a))*(1./(3*lam(a)*lam(a)+c1))*((lam(a)*lam(a)+c1)*eye(3)+lam(a)*T_flav_mat+T_sq_flav_mat);
   %disp(summ);
end

ret=summ*phi_phase;

y=ret;
end
