%Calcolo degli stati confinati
material.levelsc=solv_disp_eq1D(material.mstar,material.m2star,material.Ebulk,10*material.W,0,0); % definiti rispetto al fondo della banda di conduzione
material.levelsh=solv_disp_eq1D(material.mstarhh,material.m2starhh,material.Ebulkh,10*material.W,0,0); % definiti rispetto al fondo della banda di valenza

Ec0=material.levelsc(1);
Eh0=material.levelsh(1);

material.Egp = material.Eg + material.levelsc(1) + material.levelsh(1);

if ((length(material.levelsc)>1) & ((length(material.levelsh)>1)))
    material.Egp2 = material.Eg + material.levelsc(2) + material.levelsh(2);
else
    material.Egp2=NaN;
end

material.mrstar=(material.mstar*material.mstarhh)/(material.mstar+material.mstarhh);
material.meffel=material.mstar*const_phys.m0.data;
material.mefflac=material.mstarhh*const_phys.m0.data;



