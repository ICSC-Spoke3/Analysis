#!/bin/bash



n1=$1
n2=$2

for ((i=n1; i<=n2; i++)); do
    cat run_analysis.sh > run_analysis_$i.sh
    chmod a+x run_analysis_$i.sh
    cat Submit_clutrack.csi | sed "s/CHANGEME/$i/g" > Submit_clutrack_$i.csi 
    
done

#for i in {n1..n2} ; do cat Submit_clutrack.csi | sed "s/CHANGEME/$i/g" > Submit_clutrack_$i.csi ; done

#for i in Submit_clutrack_{${n1}..${n2}}.csi ; do 
        #echo "$i"
        #condor_submit $i -name ettore ; 

for ((i=n1;i<n2;++i));do

    echo "Submit_clutrack_$i"
    condor_submit "Submit_clutrack_$i.csi" -name ettore ;
done
