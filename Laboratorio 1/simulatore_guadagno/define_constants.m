% Physical constants used in the programs as global variables

% (la variabile globale e' una struttura in cui ciascun campo e' una costante fisica. Ogni campo 
% contiene a sua volta due sottocampi: il campo .data: valore della costante 
%                                     il campo meas_unit: unita' di misura e normalizzazione usata nei programmi)

global const_phys

const_phys.h.data=6.6261;
const_phys.h.meas_unit='*1e-34 Js Planck const';

const_phys.ht.data=6.6261/(2*pi); % h tagliato
const_phys.ht.meas_unit='*1e-34 Js Planck const';

const_phys.hteV.data=6.5821;
const_phys.hteV.meas_unit='*1e-16 eV*s Planck const';


const_phys.e.data=1.6022;
const_phys.e.meas_unit='*1e-19 C elementary charge';

const_phys.m0.data=0.91094;
const_phys.m0.meas_unit='*1e-30 Kg free electron mass';

const_phys.K.data=8.6175e-5;
const_phys.K.meas_unit='*eV/K Boltzmann constant';

const_phys.eps0.data=8.8542;
const_phys.eps0.meas_unit='*1e-12 A*s/V*m  dielectric permittivity in vacuum';

const_phys.c.data=2.9979;
const_phys.c.meas_unit='*1e8 m/s light velocity';

