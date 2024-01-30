
#!/bin/bash

#recommended code changes: g1 & g1_gas x0.01, correct e- life
#set correct saturation value PHE_MAX in NEST.hh
#S1,S2 calc modes "Waveform" in analysis.hh, and verbosity set
#to off (0 not -1) to eliminate headers from all output text

mkdir -p muonPulses/
cd muonPulses/
rm -f photon_times* track*
cd ../

for ((n=0;n<1;n++)); #number of muons to simulate
do
    #echo $n    
    for ((z=1460;z>=0;z-=1)); #top of the liq down to cath grid
    do	
	#echo $z
	#(1.122MeV-cm^2/g)(2.8811g/cm^3)(1e3keV/MeV)(0.1cm/mm)
	#=323.3 keV (MIP LET from ESTAR, LXe rho NIST)
	./execNEST 1 beta 323.3 323.3 -1 0,0,$z -1 >> muonPulses/tracker$n.txt 2> /dev/null
	mv photon_times.txt muonPulses/photon_times_$z.txt
    done
    cd muonPulses
    cat photon_times_*.txt >> photon_times$n.txt #phe not phd
    rm photon_times_*.txt
    cd ..
    
done

#current zMax above [mm] works for the LZ (or nT) detector
#replace 1st -1 in execNEST above with (uniform) field (or make z_step bigger)

#this bash loop should be the equivalent (with saturation) of
#./execNEST nMax MIP 1.122 0,0,1460 -1 0,0,0 0 2> /dev/null
