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

% X. B. 2024-10-03 Initial work

function dzdt = odefun_nucleocytoplasmictransportgrowthsensitivity1(t,z,pv,pf,phif)
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
                      28,15,1;28,16,0
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
                      28,15,1;28,1,1;28,16,0
                     ]);
% Stoichiometric matrix             
V = (N_matrix - M_matrix)';

zcyt = sum(z(1:7))+sum(z(15:16));
Vcy = pf(1)*pf(3)+pf(2)*zcyt;
NCcyto = sum(z(8:14))/zcyt;
Vn = NCcyto*Vcy;
% Constant permeability
% a = pv(14);
% Scaling with respect to nuclear surface area
% a = pv(14)*Vn^(2/3); 
% Scaling with respect to nuclear volume
a = pv(14)*Vn; 
% Experimenting with scaling
% a = pv(14)*Vn^(4/3);
% a = pv(14)*Vn^0.9;
% a = pv(14)*Vn^1.1;
% Scaling with respect to cytoplasmic protein number
% a = pv(14)*zcyt;

v = [pv(1)*(pf(3)/Vcy)*z(1)*z(2);
     pv(2)*z(3);
     pv(3)*pf(4)^(-1)*z(6);
     pv(4)*(pf(3)/Vcy)*z(7)*z(2);
     a/Vn*z(13);
     a/Vcy*z(6);
     pv(5)*(pf(3)/Vn)*z(14)*z(9);
     pv(6)*z(13);
     pv(7)*pf(4)*z(10);
     pv(8)*(pf(3)/Vn)*z(8)*z(9);
     a/Vcy*z(3);
     a/Vn*z(10);
     pv(9)*z(11);
     pv(10)*(pf(3)/Vn)*z(14)*z(12);
     a/Vcy*z(7);
     a/Vn*z(14);
     a/Vcy*z(4);
     a/Vn*z(11);
     pv(11)*(pf(3)/Vcy)*z(7)*z(5);
     pv(12)*z(4);
     a/Vn*z(8);
     a/Vcy*z(1);
     pv(13)*phif(1)*z(15);
     pv(13)*phif(2)*z(15);
     pv(13)*(1-(phif(1)+phif(2)+pv(15)+pv(16)+pv(17)))*z(15);
     pv(13)*pv(15)*z(15);
     pv(13)*pv(16)*z(15);
     pv(13)*pv(17)*z(15)
    ];
% ODE function
dzdt = V*v;
end
