syms x1 x2 real
 syms q11 q12 q13 q14 q15  real
 syms q21 q22 q23 q24 q25 real
 syms q31 q32 q33 q34 q35  real
 syms q41 q42 q43 q44 q45  real
 syms q51 q52 q53 q54 q55  real
%syms q61 q62 q63 q64 q65 q66 q67 q68 q69 real
%syms q71 q72 q73 q74 q75 q76 q77 q78 q79 real
%syms q81 q82 q83 q84 q85 q86 q87 q88 q89 real
%syms q91 q92 q93 q94 q95 q96 q97 q98 q99 real




L=[q11 q12 q13 q14 q15;% q16 q17 q18 q19
 q21 q22 q23 q24 q25;% q26 q27 q28 q29
 q31 q32 q33 q34 q35;% q36 q37 q38 q39
 q41 q42 q43 q44 q45;% q46 q47 q48 q49
 q51 q52 q53 q54 q55;% q56 q57 q58 q59
];

sigma=[x1 x2 x1^2 x1*x2 x1^3 ]'; 
%bar_sigma=[x1^2 x2^2 x1^3 x1^4 x1^5 x1^6 x1*x2 x1^2*x2 x1^3*x2]'; 
phi = [x1 x2 x1^2]';
%bar_phi=[x1^2 x1*x2 x2^2 x1^3 x1^2*x2 x1^4]';



 syms p11 p12 p13   real
 syms p21 p22 p23  real
 syms p31 p32 p33 real

 
 P=[p11 p12 p13  
    p21 p22 p23  
    p31 p32 p33 ];

% V=phi'*P*phi
% 
% Vx1=collect(diff(V,x1),[x1,x2])
% Vx2=collect(diff(V,x2),[x1,x2])


 %P=magic(3)
 %sigma=[x1 x2 x1^2 x1*x2 x1^3 ]'; 
W=[2*P(1,1) P(1,2)+P(2,1) 3*P(1,3)+3*P(3,1) 2*P(2,3)+2*P(3,2) 4*P(3,3);
  P(1,2)+P(2,1) 2*P(2,2) P(2,3)+P(3,2) 0 0];     
% ]

F=[0 -1 -3/2 0 -1/2;
   0 0   0   0 0];%
G=[0;1];
Q=diag([.1,0,0,0, 0]);
R=1;
K=[6 -3 0 0 0];%

collect(sigma'*L*sigma,[x1,x2])

%diff(phi'*P*phi,x2)-W(2,:)*sigma
