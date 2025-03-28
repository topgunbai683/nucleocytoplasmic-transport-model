% ODE function for law of mass action approximation of a simplified
% model of Ran-mediated nucleocytoplasmic transport
% coupled with translation and growth.
% Based on: Wang et al. Thermodynamic Paradigm for Solution Demixing 
% Inspired by Nuclear Transport in Living Cells,
% and our previous work.
% Cytoplasmic volume Vcy is set to be an affine function of the total number 
% of cytoplasmic proteins.
% Nuclear volume Vn is variable and is determined by the osmotic balance
% at the nuclear envelope.
% dzdt is calculated as the matrix multiplication of the stoichiometric matrix
% and the propensity function vector.
% Assume Cytoplasmic RanGDP can diffuse through NPCs with a much lower
% diffusion rate than others.

% X. B. 2025-02-20 Initial work

function dzdt = odefun_nucleocytoplasmictransportgrowthran1(t,z,p)
% Reaction matrices
M_matrix = spconvert([1,1,1;1,2,1;
                      2,3,1;
                      3,6,1;
                      4,2,1;4,7,1;
                      5,13,1;
                      6,6,1;
                      7,9,1;7,14,1;
                      8,13,1;
                      9,10,1;
                      10,8,1;10,9,1;
                      11,3,1;
                      12,10,1;
                      13,11,1;
                      14,12,1;14,14,1;
                      15,7,1;
                      16,14,1;
                      17,4,1;
                      18,11,1;
                      19,5,1;19,7,1;
                      20,4,1;
                      21,8,1;
                      22,1,1;
                      23,15,1;
                      24,15,1;
                      25,15,1;
                      26,15,1;
                      27,15,1;
                      28,15,1;
                      29,2,1;
                      30,17,1;
                      31,17,1;
                      32,9,1;32,17,0
                     ]);
N_matrix = spconvert([1,3,1;
                      2,1,1;2,2,1;
                      3,2,1;3,7,1;
                      4,6,1;
                      5,6,1;
                      6,13,1;
                      7,13,1;
                      8,9,1;8,14,1;
                      9,8,1;9,9,1;
                      10,10,1;
                      11,10,1;
                      12,3,1;
                      13,12,1;13,14,1;
                      14,11,1;
                      15,14,1;
                      16,7,1;
                      17,11,1;
                      18,4,1;
                      19,4,1;
                      20,5,1;20,7,1;
                      21,1,1;
                      22,8,1;
                      23,15,2;
                      24,15,1;24,5,1;
                      25,15,1;25,16,1;
                      26,15,1;26,7,1;
                      27,15,1;27,2,1;
                      28,15,1;28,1,1;
                      29,17,1;
                      30,2,1;
                      31,9,1;
                      32,17,1
                     ]);
% Stoichiometric matrix             
V = (N_matrix - M_matrix)';

zcyt = sum(z(1:7))+sum(z(15:16));
Vcy = p.V0+p.Ccy*zcyt;
NCcyto = (sum(z(8:14))+z(17))/zcyt;
Vn = NCcyto*Vcy;
% Constant permeability
% a = p.a0;
% Scaling with respect to nuclear surface area
% a = p.Cnp*Vn^(2/3); 
% Scaling with respect to nuclear volume
a = p.Cnp*Vn; 
% Experimenting with scaling
% a = p.Cnp*Vn^(4/3);
% a = p.Cnp*Vn^0.9;
% a = p.Cnp*Vn^1.1;
% Scaling with respect to cytoplasmic protein number
% a = p.Cnp*zcyt;

v = [p.k1f*(p.Vref/Vcy)*z(1)*z(2);
     p.k1b*z(3);
     p.k2f*z(6); % Constant
     p.k2b*(p.Vref/Vcy)*z(7)*z(2);
     a/Vn*z(13);
     a/Vcy*z(6);
     p.k4f*(p.Vref/Vn)*z(14)*z(9);
     p.k4b*z(13);
     p.k5f*z(10); % Constant
     p.k5b*(p.Vref/Vn)*z(8)*z(9);
     a/Vcy*z(3);
     a/Vn*z(10);
     p.k7f*z(11);
     p.k7b*(p.Vref/Vn)*z(14)*z(12);
     a/Vcy*z(7);
     a/Vn*z(14);
     a/Vcy*z(4);
     a/Vn*z(11);
     p.k10f*(p.Vref/Vcy)*z(7)*z(5);
     p.k10b*z(4);
     a/Vn*z(8);
     a/Vcy*z(1);
     p.k12f*z(15);
     p.k13f*z(15);
     p.k14f*z(15);
     p.k15f*z(15);
     p.k16f*z(15);
     p.k17f*z(15);
     p.Cran*a/Vcy*z(2);
     p.Cran*a/Vn*z(17);
     p.k19f*z(17);
     p.k19b*z(9)
    ];
% ODE function
dzdt = V*v;
end
