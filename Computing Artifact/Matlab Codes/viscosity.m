%---------------------------------------------------------------------%
%This code computes the value of nu - constant term.
%Written by Olayemi Adeyemi on 11/2022
%           Computing Artifact
%           Computing PhD 
%           Boise State University
%---------------------------------------------------------------------%


%To obtain the value of nu - viscosity term
    sigma_t = 1;
    c_mu = 0.09;
    rho = 1.204;
    u = 2;
    L = 1; %mountain height;
    Re = 10e6;
    I = 0.16*(Re^(-0.125));
    e = 1.5*((u*I)^2);
    l = 0.07*L;
    epsilon = ((c_mu)^(0.75))*(e^1.5)*(l^(-1));
    mu_t = c_mu*rho*(e^2 / epsilon);
    nu = mu_t / sigma_t; 
    nu;