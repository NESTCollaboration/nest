
#!/bin/bash

for e in 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.5 0.52 0.54 0.56 0.58 0.6 0.62 0.64 0.66 0.68 0.7 0.72 0.74 0.76 0.78 0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 1 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.75 1.8 1.85 1.9 1.95 2 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5 5.2 5.4 5.6 5.8 6 6.2 6.4 6.6 6.8 7 7.2 7.4 7.6 7.8 8 8.2 8.4 8.6 8.8 9 9.2 9.4 9.6 9.8 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 50.5 51 51.5 52 52.5 53 53.5 54 54.5 55 55.5 56 56.5 57 57.5 58 58.5 59 59.5 60 60.5 61 61.5 62 62.5 63 63.5 64 64.5 65 65.5 66 66.5 67 67.5 68 68.5 69 69.5 70 71 72 73 74 75 76 77 78 79 80 82 84 86 88 90 95 100
#add lower Es even as low as 0.1 (keV) if doing ER, and reduce high end to ~20-40

do
    #Uncomment the echo statement and comment the execution line out to get the energies printed out, to e.g. add to a text output file after the fact for better organization
    #echo $e
    /Users/szydagis/Desktop/buildNEST/execNEST 1e5 NR $e $e 192 -1 0 > /dev/null
    #Don't forget to change the above to your executable's absolute path, and change NR to beta for doing ER instead. Change 192 to your detector's drift E-field (in V/cm)
    #1e5: note 1e6 is also a reasonable number of events. Make sure to use the "2>" operator to redirect the terminal screen output of this script to a text file to save it.

done

#Lastly, copy/paste results from Excel into txt in terminal, using nest/examples/LUXRun03_nrEff_Simulated.txt as example
