close all ;clear ;
%% Given data

rho= 1.225; % Density of air
Nb= 4; % Number of blades
R= 6; % Blade radius
c= 0.4; % Blade chord
Cd0= 0.01; % Profile Drag coefficient
Cla= 5.73; % Lift curve slope
omega= 10*pi; % Rotor angular rate
R_ct= 0.2*R; % Root cut-out 
CT_h= 0.006; % Hover CT 
n= 10; % number of radial stations/segments

%% calculated data

r_ct= 0.2; % Non-dimensional root cutout
sigma= Nb*c*(R-R_ct)/(pi*R^2); % Rotor Solidity
lambda_h= sqrt(CT_h/2); % induced inflow for hover
dy= (1-0.2)/n; % Non-dimensional length of each segment
theta_0= (6*CT_h/sigma/Cla)+(3*sqrt(CT_h/2))/2; % reference pitch angle from BET&MT
theta_75=theta_0; % reference pitch angle considered at 0.75R
theta_tip= (4*CT_h/sigma/Cla) + sqrt(CT_h/2); % Ideal Twist solution for pitch angle
c_tip= 0.5*c; % Ideal Taper assumption
sigma_tip= Nb*c_tip*(R-R_ct)/(pi*R^2); % rotor solidity for ideal taper
% Gaussian Quadrature Coefficients

weight= [0.171324492 0.360761573 0.467913935 0.467913935 0.360761573 0.171324492];
node=[-0.9324695142031520278123,-0.661209386466264513661,-0.2386191860831969086305,....
    0.238619186083196908631,0.661209386466264513661,0.9324695142031520278123];


%% Numerical Values Initialization

k=1; % array initializer for mean values
CT_i= zeros(1,n);
CT_ii= zeros(1,n);
CT_iii= zeros(1,n);
CT_iv= zeros(1,n);

CPi_i= zeros(1,n);
CPi_ii= zeros(1,n);
CPi_iii= zeros(1,n);
CPi_iv= zeros(1,n);

f_ct_i= zeros(1,6);
f_ct_ii= zeros(1,6);
f_ct_iii= zeros(1,6);
f_ct_iv= zeros(1,6);

f_cp_i= zeros(1,6);
f_cp_ii= zeros(1,6);
f_cp_iii= zeros(1,6);
f_cp_iv= zeros(1,6);
CP0_iv= zeros(1,n);

lambda_collect_i= 0;
lambda_collect_ii= 0;
lambda_collect_iii= 0;
lambda_collect_iv= 0;

for i=1:1:n % span devided into n number of radial stations/segments
    
    %%
    % _% Gaussian Quadrature for each segment_ 
    a= 0.2+dy*(i-1); % integration lower limit
    b= 0.2+dy*i; % integration upper limit
    t= ((b-a)/2)*node + ((b+a)/2); % change of variables: new locations are t(i) from x(i)
    theta_tw_ii=0; % for question 1-(ii) theta twist rate is zero
    theta_ii= theta_75 + theta_tw_ii*(t-0.75); % theta for zero twist rate
    theta_tw_iii=-15*pi/180; % for question 1-(iii) theta twist rate is -15deg
    theta_iii= theta_75 + theta_tw_iii*(t-0.75); % theta for -15deg twist rate
    sigma_iv= Nb*(c_tip./t)*(R-R_ct)/(pi*R^2); % for question 1-(iv) ideal taper
    for j=1:1:6
        %% C1. Ideal Twist: theta_tip

        theta_tip_i= theta_tip; % theta_tip calculation
        lambda_t_i= (sigma*Cla/16)*(sqrt(1+(32*theta_tip_i/sigma/Cla))-1); % lambda for each loaction t
        lambda_collect_i=  lambda_collect_i + lambda_t_i; % collecting all lambda over one gaussian span for mean
        dCT_dr_i(j)= 0.5*sigma*Cla*(theta_tip_i-lambda_t_i)*t(j); % CT distribution from lambda
        f_ct_i(j)= dCT_dr_i(j); % gaussian function for CT
        CT_i(i)= CT_i(i) + weight(j)*f_ct_i(j)*0.5*(b-a); % numerically performing integration CT over t
        f_cp_i(j)= 0.5*sigma*Cla*(theta_tip_i-lambda_t_i)*t(j)*lambda_t_i; % gaussian integration for CPi, induced power
        CPi_i(i)= CPi_i(i) + weight(j)*f_cp_i(j)*0.5*(b-a); % numerically performing integrating CPi over t
        
        %% C2. No twist: theta_tw= 0deg

        theta_tip_ii= theta_ii(j)*t(j);
        lambda_t_ii= (sigma*Cla/16)*(sqrt(1+(32*theta_tip_i/sigma/Cla))-1);
        lambda_collect_ii=  lambda_collect_ii + lambda_t_ii;
        dCT_dr_ii(j)= 0.5*sigma*Cla*(theta_tip_ii-lambda_t_ii)*t(j);
        f_ct_ii(j)= dCT_dr_ii(j);
        CT_ii(i)= CT_ii(i) + weight(j)*f_ct_ii(j)*0.5*(b-a);
        f_cp_ii(j)= 0.5*sigma*Cla*(theta_tip_ii-lambda_t_ii)*t(j)*lambda_t_ii;
        CPi_ii(i)= CPi_ii(i) + weight(j)*f_cp_ii(j)*0.5*(b-a);
        
        %% C3. Linear Twist: theta_tw= -15deg

        theta_tip_iii= theta_iii(j)*t(j);
        lambda_t_iii= (sigma*Cla/16)*((1+(32*theta_tip_iii/sigma/Cla))^0.5-1);
        lambda_collect_iii=  lambda_collect_iii + lambda_t_iii;
        dCT_dr_iii(j)= 0.5*sigma*Cla*(theta_tip_iii-lambda_t_iii)*t(j);
        f_ct_iii(j)= dCT_dr_iii(j);
        CT_iii(i)= CT_iii(i) + weight(j)*f_ct_iii(j)*0.5*(b-a);
        f_cp_iii(j)= 0.5*sigma*Cla*(theta_tip_iii-lambda_t_iii)*t(j)*lambda_t_iii;
        CPi_iii(i)= CPi_iii(i) + weight(j)*f_cp_iii(j)*0.5*(b-a);
        
        %% C4. Ideal Twist: theta_tip and Ideal Taper: c_tip= 0.5*c

        theta_tip_iv= theta_tip;
        lambda_t_iv= (sigma_iv(j)*Cla/16)*((1+(32*theta_tip_iv/sigma_iv(j)/Cla))^0.5-1);
        lambda_collect_iv=  lambda_collect_iv + lambda_t_iv;
        dCT_dr_iv(j)= 0.5*sigma_iv(j)*Cla*(theta_tip_iv-lambda_t_iv)*t(j);
        f_ct_iv(j)= dCT_dr_iv(j);
        CT_iv(i)= CT_iv(i) + weight(j)*f_ct_iv(j)*0.5*(b-a);
        f_cp_iv(j)= 0.5*sigma_iv(j)*Cla*(theta_tip_iv-lambda_t_iv)*t(j)*lambda_t_iv;
        CPi_iv(i)= CPi_iv(i) + weight(j)*f_cp_iv(j)*0.5*(b-a);
        f_cp0_iv(j)= 0.5*sigma_iv(j)*Cd0*((t(j))^3);
        CP0_iv(i)= CP0_iv(i) + weight(j)*f_cp0_iv(j)*0.5*(b-a);
    end
    
    loc(k)= t(3); % for each gaussian span it is a middle location: for plotting the mean values over the entire span
    %% mean values for each element Case-1

    lambda_mean_i(k)= lambda_collect_i/6;            % induced inflow
    lambda_collect_i=0;
    CT_i_mean(k)= mean(dCT_dr_i);                    %thrust distribution
    CPi_i_mean(k)= lambda_mean_i(k)*CT_i_mean(k);    % induced power
    CP0_i_mean(k)= sigma*Cd0*((t(3))^3)/2;           % profile power
    CP_i_mean(k)= CPi_i_mean(k) + CP0_i_mean(k);     % total power
    %% mean values for each element Case-2

    lambda_mean_ii(k)= lambda_collect_ii/6;
    lambda_collect_ii=0;
    CT_ii_mean(k)= mean(dCT_dr_ii);
    CPi_ii_mean(k)= lambda_mean_ii(k)*CT_ii_mean(k);
    CP0_ii_mean(k)= sigma*Cd0*((t(3))^3)/2;
    CP_ii_mean(k)= CPi_ii_mean(k) + CP0_ii_mean(k);
    %% mean values for each element Case-3

    lambda_mean_iii(k)= lambda_collect_iii/6;
    lambda_collect_iii=0;
    CT_iii_mean(k)= mean(dCT_dr_iii);
    CPi_iii_mean(k)= lambda_mean_iii(k)*CT_iii_mean(k);
    CP0_iii_mean(k)= sigma*Cd0*((t(3))^3)/2;
    CP_iii_mean(k)= CPi_iii_mean(k) + CP0_iii_mean(k);
    %% mean values for each element Case-4

    lambda_mean_iv(k)= lambda_collect_iv/6;
    lambda_collect_iv=0;
    CT_iv_mean(k)= mean(dCT_dr_iv);
    CPi_iv_mean(k)= lambda_mean_iv(k)*CT_iv_mean(k);
    CP0_iv_mean(k)= sigma_iv(3)*Cd0*((t(3))^3)/2;
    CP_iv_mean(k)= CPi_iv_mean(k) + CP0_iv_mean(k);
    
     k= k+1; % Incrementing the location
    
end


%% Results and Plots

%% plot 1  'Pitch Angle vs. r'
figure('Name','Pitch Angle vs. r')
r= r_ct:r_ct/10:1;

% Ideal Twist
CT_i_total= sum(CT_i); % total CT over the entire disc
CP_i_total= sum(CPi_i)+(sigma*Cd0/8); % total CP over the entire disc
lambda_i= (CT_i_total/2)^0.5; % recalculating induced inflow from final CT
theta_75_i= ((6*CT_i_total/sigma/Cla)+(3*lambda_i)/2);
theta_tip_IdealTwist= (4*CT_i_total/sigma/Cla) + lambda_i; % Theta_tip closed form solution
theta_result_i= (theta_tip_IdealTwist./r)*(180/pi); % theta vs r plot


% No twist: theta_tw= 0deg
CT_ii_total= sum(CT_ii);
CP_ii_total= sum(CPi_ii)+(sigma*Cd0/8);
lambda_ii= (CT_ii_total/2)^0.5;
theta_result_ii= ((6*CT_ii_total/sigma/Cla)+(3*lambda_ii)/2)*(180/pi).*ones(1,length(r));
theta_tip_NoTwist= theta_result_ii.*r;


% Linear Twist: theta_tw= -15deg
CT_iii_total= sum(CT_iii);
CP_iii_total= sum(CPi_iii)+(sigma*Cd0/8);
lambda_iii= (CT_iii_total/2)^0.5;
theta_0_iii= ((6*CT_iii_total/sigma/Cla)+(3*lambda_iii)/2);
theta_result_iii= (theta_0_iii + theta_tw_iii*(r-0.75))*(180/pi);
theta_tip_Twist= theta_result_iii.*r;


% Ideal Twist: theta_tip and Ideal Taper: c_tip= 0.5*c
CT_iv_total= sum(CT_iv);
CP_iv_total= sum(CPi_iv)+sum(CP0_iv);
lambda_iv= (CT_iv_total/2)^0.5;
sigma_75= sigma_tip/0.75;
theta_75_iv= ((6*CT_iv_total/sigma_75/Cla)+(3*lambda_iv)/2);
theta_result_iv= (((4*CT_iv_total/sigma_tip/Cla)+lambda_iv)./r)*(180/pi);
theta_tip_IdealTaper= (4*CT_iv_total/sigma_tip/Cla)+lambda_iv;
plot(r,theta_result_i,'b--',r,theta_result_ii,'r',r,theta_result_iii,'g.-',r,theta_result_iv,'y-*')
title('Variation of Pitch angle vs. Non-dimentional radial location Curves')
xlabel('Non-dimentional radial location, r')
ylabel('Pitch angle, Theta')
xlim([0.19,1])
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
%% plot 2 'Numerical Solution vs. Analytical Solution'

figure('Name','Numerical Solution vs. Analytical Solution','NumberTitle','off')
plot(r,theta_result_i,'y--',r,(theta_tip./r),'g.-') % Numerical Solution, Analytical formula
title('Comparison of Numerical solution with closed form exact formula for ideal twist')
xlabel('Non-dimentional radial location, r')
ylabel('Pitch angle, Theta')
xlim([0.18,1])
legend('Numerical Solution','Closed form exact solution')

%% plot 3  Angle of attack vs non dim. radial location

figure('Name','Variation of AOA vs. r','NumberTitle','off')
% Ideal Twist: theta_tip
alpha_i= theta_result_i-(lambda_i./r); % Calculating AOA for all the r locations
% No twist: theta_tw= 0deg
alpha_ii= theta_result_ii-(lambda_ii./r);
% Linear Twist: theta_tw= -15deg
alpha_iii= theta_result_iii-(lambda_iii./r);
% Ideal Twist: theta_tip and Ideal Taper: c_tip= 0.5*c
alpha_iv= theta_result_iv-(lambda_iv./r);
plot(r,alpha_i,'b--',r,alpha_ii,'r.',r,alpha_iii,'g',r,alpha_iv,'y-*') 
title('AOA vs r','FontSize',15)
xlabel('Non-dimentional radial location, r')
ylabel('AOA, (alpha)')
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
xlim([0.18,1])


%% plot 4 Inflow vs r

figure('Name','Lambda vs. r','NumberTitle','off')
plot(loc,lambda_mean_i,'b-|',loc,lambda_mean_ii,'r.-',loc,lambda_mean_iii,'g-o',loc,lambda_mean_iv,'y-*') % plotting the induced inflow for all r locations
title('Variation of Induced Inflow vs. Non-dimentional radial location Curves')
xlabel('Non-dimentional radial location, r')
ylabel('Induced Inflow, Lambda')
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
xlim([0.18, 1])

%% plot 5 dCT/dr vs r

figure('Name','dCT/dr vs. r','NumberTitle','off')
plot(loc,CT_i_mean,'b--',loc,CT_ii_mean,'r.-',loc,CT_iii_mean,'g-o',loc,CT_iv_mean,'y-*') % plotting the CT distribution for all r locations
title('Variation of thrust distribution vs. Non-dimentional radial location Curves')
xlabel('Non-dimentional radial location, r')
ylabel('thrust distribution, dCT/dr')
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
xlim([0.18, 1])

%% plot 6 dCQi/dr vs. r

figure('Name','dCQi/dr vs. r','NumberTitle','off')
plot(loc,CPi_i_mean,'b--',loc,CPi_ii_mean,'r.-',loc,CPi_iii_mean,'g-o',loc,CPi_iv_mean,'y-*') % plotting the induced power distn. for all r locations
title('Variation of induced torque distribution vs. Non-dimentional radial location Curves')
xlabel('Non-dimentional radial location, r')
ylabel('induced torque distribution, dCQi/dr')
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
xlim([0.18, 1])

%% plot 7 dCQ0/dr vs. r

figure('Name','dCQ0/dr vs. r','NumberTitle','off')
plot(loc,CP0_i_mean,'b--',loc,CP0_ii_mean,'r.-',loc,CP0_iii_mean,'g-o',loc,CP0_iv_mean,'y-*') % plotting the profile power distn. for all r locations
title('Variation of profile torque distribution vs. Non-dimentional radial location Curves')
xlabel('Non-dimentional radial location, r')
ylabel('profile torque distribution, dCQ0/dr')
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
xlim([0.18,1])

%% plot 8 dCQ/dr vs. r

figure('Name','dCQ/dr vs. r')
plot(loc,CP_i_mean,'b--',loc,CP_ii_mean,'r.-',loc,CP_iii_mean,'g-o',loc,CP_iv_mean,'y-*') % plotting the total power distn. for all r locations
title('Total torque distribution vs. Non-dimentional radial location')
xlabel('Non-dimentional radial location, r')
ylabel('total torque distribution, dCQ/dr')
legend('Ideal Twist: theta_tip','No twist: theta_tw= 0deg'...
    ,'Linear Twist: theta_tw= -15deg','Ideal Twist: theta_tip and Ideal Taper')
xlim([0.18 1])
