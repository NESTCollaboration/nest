
#!/bin/bash

for a in -0.1 0.00 0.10
do
for b in 0.40 0.50 0.60
do
for c in 0.04 0.05 0.06
do
for d in -0.5 -0.4 -0.3
do
for e in 1.09 1.10 1.11
do
for f in 0.96 0.97 0.98
do
for g in 6e-2 7e-2 8e-2
do
    
    echo $a $b $c $d $e $f $g
    #/home/mszydagis/build/testNEST $a $b $c $d $e $f $g > /home/mszydagis/nest_output/nest_output.ER 2> /dev/null
    #/home/mszydagis/build/rootNEST /home/mszydagis/nest_output/nest_output.ER /home/mszydagis/build/Run4_C14_Time_bin_5_105_tdrift_170_us.txt 2> /dev/null
    #rm /home/mszydagis/nest_output/nest_output.ER
    
done
done
done
done
done
done
done
