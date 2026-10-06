Network of autoscopic hallucinations elicited by intracerebral stimulations of periventricular nodular heterotopia: an SEEG study 
*Flavius-Ionut Bratua, *Irina Oanea, Andrei Barborica, Cristian Donos, Constantin Pistol, Andrei Daneasa, Camelia Lentoiu, Ioana Mindruta

This folder contains the data and scripts required to reproduce the results presented in the paper.

Folders:
	- MRI - anonymized T1 MRI of the patient
	- Fibers - contains the source (*.src) and fiber (*.fib) files, as well as ROIs and tracts used by DSI Studio to generate figures 1 and 2
	- Matlab - Matlab scripts for calculating the three distinct connectivities.
		- seeg_fibers.m - used for calculating anatomical connectivity
		- pvnh_h2.m - used for calculating the h2 and SPES connectivity
	- Data - data files
		- MRI_defaced.mgz defaced and anonymized MRI
		- HFS - 43 Hz stimulations in AnyWave (*.ades) format.
		- SPES - single-pulse electrical stimulation data in AnyWave (*.ades) format.

The ADES format specification is available at https://meg.univ-amu.fr/wiki/AnyWave:ADES 