
#!/bin/bash
#This code doesn't actually perform any model fitting, however,
#it provides an example for using rootNEST if one wishes to use testNEST + rootNEST
#to fit a data target to find optimal model/detector parameters to match data

#See lines 70-140 of src/testNEST.cpp for more info on what each of these input arguments is doing

#Note that rootNEST must be in FIT mode (mode=1) to compare with data

#First declare ranges for each of the free parameters/
#Here, in several for-loops, using 7 parameters with set ranges
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
    
    #below: run testNEST with the select parameters, and save the output
    #/path/to/testNEST $a $b $c $d $e $f $g > nest_output.ER 2> /dev/null
    
    #next, use rootNEST in FIT mode to get goodness-of-fit statistic with a data file (See README on rootNEST)
    #/path/to/rootNEST nest_output.ER /path/to/fittingTarget_dataFile.dat 2> /dev/null
    
    #then remove the output file and move on to the next parameter set
    #rm nest_output.ER
    
done
done
done
done
done
done
done
