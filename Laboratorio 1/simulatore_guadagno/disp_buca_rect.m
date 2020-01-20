function disp_eq=disp_buca_rect(E,m1,m2,DE,w,n_soluz)
global const_phys

% structure, reference levels, input parameters:
%
%	POTENZIAL:
%
%  --------------                 ----------------  energy level=DE     
%						|						  |							  
%               | effective mass: |   effective mass:  
%						|	m1					  |   m2 
%						|                 |
%						------------------ energy level=0
% 
% inputs: 
% E: energy level of the confined state
% m1, m2: effective masses in the two regions (well and barrier) normalized 
% respect to the free electron  mass
% DE: depth of the well in eV
% w: width of the well in Amstrong
% n_soluz: solution number


norm_const=0.38203508582328; % normalization constant 
disp_eq=(sqrt(2*E*m1)*w/2/const_phys.ht.data*norm_const-(n_soluz)*pi/2)-atan(sqrt(m1/m2*(DE-E)/E));
   
