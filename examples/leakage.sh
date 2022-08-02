
#!/bin/bash

for f in 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400 450 500 550 600 700 800 900 1000 1200 1400 1600 1800 2000 2500 3000 3500 4000 4500 5000 5500 6000

do

    ./execNEST 5e6 beta 0 25 $f -1 0 > NESTOutputER.new
    ./execNEST 2e6 D-D 0 74 $f 200,0,-1 0 > NESTOutputNR.new
    examples/rootNEST NESTOutputER.new NESTOutputNR.new

done
