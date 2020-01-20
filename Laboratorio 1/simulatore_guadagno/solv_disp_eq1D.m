function Elevel=solv_disp_eq1D(m1,m2,DE,w,flag_max,nmax)
global const_phys

% solution of the dispersion equation for a rectangular potential barrier using the function fzero of MATLAB
% for comments on input paramters see disp_buca_rect.m
% flag_max: =0: calculates all the energy levels
% flag_max: =1: calculates the energy levels between 0 and nmax


% returns the vector Elevel containing the energy levels of all the confined states


norm_const1=6.85162365268953;
norm_const2=0.38203508582328;

% maximum order of the solution:
% the order of the solutions are 0,1,2,...
n_soluz_max=floor(norm_const2*w*sqrt(2*m1*DE)/(pi*const_phys.ht.data));


if (flag_max==1)&(n_soluz_max>nmax)
   n_soluz_max=nmax;
end



options=optimset('fzero'); % (matlab command) set options paramters for the fzero function
options=optimset(options,'display','off');

for n_soluz=0:n_soluz_max
	% [Emin Emax] energy range where dispersion eq. changes sign
   Emax=norm_const1*1/(2*m1)*(const_phys.ht.data*(n_soluz+1)*pi/w)^2;
	if (Emax>DE)
   	Emax=DE;
	end
	Emin=norm_const1*1/(2*m1)*(const_phys.ht.data*(n_soluz)*pi/w)^2;
	if (Emin==0)
   	Emin=1e-10;
	end
   %[En,fval,exitflag,output]=fzero('disp_buca_rect',[Emin Emax],options,m1,m2,DE,w,n_soluz);
   En=fzero('disp_buca_rect',[Emin Emax],options,m1,m2,DE,w,n_soluz);
   Elevel(n_soluz+1)=En;
end

