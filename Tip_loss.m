% Problem 2
clear all
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
n= 100; % number of radial stations/segments

%% calculated data

r_ct= 0.2; % Non-dimensional root cutout
sigma= Nb*c*(R-R_ct)/(pi*R^2); % Rotor Solidity
lambda_h= sqrt(CT_h/2); % induced inflow for hover
dy= (1-0.2)/n; % Non-dimensional length of each segment
theta_0= (6*CT_h/sigma/Cla)+(3*sqrt(CT_h/2))/2; % reference pitch angle from BET&MT
error= 10^-6; % error acceptance value
% Gaussian Quadrature Coefficients

weight= [0.171324492 0.360761573 0.467913935 0.467913935 0.360761573 0.171324492];
node=[-0.9324695142031520278123,-0.661209386466264513661,-0.2386191860831969086305,....
    0.238619186083196908631,0.661209386466264513661,0.9324695142031520278123];


%% Numerical Solution

% Initialization of values
k=1; % array initializer for plots
CT(n)= 0;
CPi(n)= 0;
f_ct(6)= 0;
f_cp(6)= 0;

CT_TL(n)= 0; % while considering tip loss effect 
CPi_TL(n)= 0;
f_ct_TL(6)= 0;
f_cp_TL(6)= 0;

lambda_collect= 0;              % inflow analytical
lambda_collect_TL= 0;           % inflow using Prandtl's function

for i=1:1:n % span devided into n number of radial stations/segments
    %% Gaussian Quadrature for each segment
    
    a= 0.2+dy*(i-1); % integration lower limit
    b= 0.2+dy*i; % integration upper limit
    t= ((b-a)/2)*node + ((b+a)/2); % change of variables: new locations are t(i) from x(i)
    theta_tw=0; % no twist
    theta= theta_0 + theta_tw*(t-0.75); % final theta for initial guess
    for j=1:1:6
        % No twist: theta_tw= 0deg
        theta_tip= theta(j)*t(j); % calculating theta_tip
        lambda_t= (sigma*Cla/16)*((1+(32*theta(j)*t(j)/sigma/Cla))^0.5-1); % inflow varying with r
        e=lambda_h; % error initialization
        lambda_TL_t= lambda_t; % initially starting with F=1;
        while e>error
                f= (Nb/2)*(1-t(j))/lambda_TL_t;
                F= (2/pi)*acos(exp(-f)); % tip loss coefficients
                lambda_clc= (sigma*Cla/16/F)*((1+(32*F*theta_tip/sigma/Cla))^0.5-1); % prandtl's tip loss function
                e= abs(lambda_TL_t-lambda_clc); % error calculation
                lambda_TL_t= lambda_clc; %  lambda for next iteration
        end
        lambda_collect=  lambda_collect + lambda_t; % collecting lambda for plotting
        lambda_collect_TL=  lambda_collect_TL + lambda_TL_t;
        dCT_dr(j)= 0.5*sigma*Cla*(theta_tip-lambda_t)*t(j); % CT distribution
        dCT_dr_TL(j)= 0.5*sigma*Cla*(theta_tip-lambda_TL_t)*t(j);
        
        f_ct(j)= dCT_dr(j);
        f_ct_TL(j)= dCT_dr_TL(j);
        
        CT(i)= CT(i) + weight(j)*f_ct(j)*0.5*(b-a); % performing numerical integration for CT
        CT_TL(i)= CT_TL(i) + weight(j)*f_ct_TL(j)*0.5*(b-a);
        
        f_cp(j)= (0.5*sigma*Cla*(theta_tip-lambda_t)*t(j)*lambda_t)+(0.5*sigma*Cd0*(t(j)^3));
        f_cp_TL(j)= (0.5*sigma*Cla*(theta_tip-lambda_TL_t)*t(j)*lambda_TL_t)+(0.5*sigma*Cd0*(t(j)^3));
        
        CPi(i)= CPi(i) + weight(j)*f_cp(j)*0.5*(b-a); % performing numerical integration for CPi
        CPi_TL(i)= CPi_TL(i) + weight(j)*f_cp_TL(j)*0.5*(b-a);
    end
    
    loc(k)= t(3); % locations for plot
    lambda_mean(k)= lambda_collect/6;        % Inflow analytical  
    lambda_mean_TL(k)= lambda_collect_TL/6;  % Inflow using Prandtl's function
    
    lambda_collect=0; % collect values to zero for next iteration
    lambda_collect_TL=0;
    
    % getting the plot variables for CT, CPi for both tiploss and no
    % tiploss
    CT_mean(k)= mean(dCT_dr);
    CT_TL_mean(k)= mean(dCT_dr_TL);
    
    CP_mean(k)=CPi(i);
    CP_TL_mean(k)=CPi_TL(i);
    k= k+1; % Incrementing the location
end
%% Results and Plots

%% plot 1: pitch  vs. r
figure('Name','Pitch refrence Angle theeta vs. r')
r= linspace(r_ct,1,100); 

% No twist: theta_tw= 0deg
CT_total= sum(CT);
CT_TL_total= sum(CT_TL); 

lambda= (CT_total/2)^0.5;
lambda_TL= (CT_TL_total/2)^0.5;

theta75_result= ((6*CT_total/sigma/Cla)+(3*lambda)/2).*ones(1,length(r));
theta75_result_TL= ((6*CT_TL_total/sigma/Cla)+(3*lambda_TL)/2).*ones(1,length(r));

plot(r,theta75_result,'r-|',r,theta75_result_TL,'k-','LineWidth',1.5);
grid on
title('Pitch refrence ( Theeta_7_5 ) angle vs. r')
xlabel('Non-dimentional radial location, r')
ylabel('Pitch angle, Theta')
ylim([0.15 0.20])
legend('No twist No Tip-loss','No Twist With Tip-loss')

%% plot 2: lambda vs. r
figure('Name','Lambda vs. r','NumberTitle','off')
plot(loc,lambda_mean,'k-|',loc,lambda_mean_TL,'r-','LineWidth',1.5)
title('Induced Inflow vs. r')
xlabel('Non-dimentional radial location, r')
ylabel('Induced Inflow, Lambda')
legend('No twist No Tip-loss','No Twist With Tip-loss')
grid on

%% Plot 3: CT distn. vs. r 
figure('Name','dCT/dr vs. r')
plot(loc,CT_mean,'k-|',loc,CT_TL_mean,'r-','LineWidth',1.5)
title('Thrust distribution vs. r')
xlabel('Non-dimentional radial location, r')
ylabel('thrust distribution, CT')
legend('No twist No Tip-loss','No Twist With Tip-loss')
grid on

%% Plot 4: CP distribution Vs. r
figure('Name','dCQ/dr vs. r')
plot(loc,CP_mean,'k-|',loc,CP_TL_mean,'r-','LineWidth',1.5)
title('Torque distribution vs. r')
xlabel('Non-dimentional radial location, r')
ylabel('total torque distribution, CQ')
legend('No twist No Tip-loss','No Twist With Tip-loss')
grid on
