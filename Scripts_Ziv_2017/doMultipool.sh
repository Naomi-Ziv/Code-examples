#!/bin/sh

#all folders have same file names (chr01.txt, chr02.txt...)
cd ~/multipool-master/C1S2
#loop on chromosomes
for i in *txt
do
    #loop on different recombination fraction
    for j in 1000 2500
    do
	#loop on different population size 
	for k in 200 1000
	do

	    #C1
	    echo
	    echo C1_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/C1.${i}.${j}.${k} -np ~/multipool-master/C1S2/${i} ~/multipool-master/C1S6/${i}

	    #C2
	    echo
	    echo C2_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/C2.${i}.${j}.${k} -np ~/multipool-master/C2S2/${i} ~/multipool-master/C2S6/${i}

	    #C3
	    echo
	    echo C3_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/C3.${i}.${j}.${k} -np ~/multipool-master/C3S2/${i} ~/multipool-master/C3S6/${i}

	    #C4
	    echo
	    echo C4_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/C4.${i}.${j}.${k} -np ~/multipool-master/C4S1/${i} ~/multipool-master/C4S6/${i}

	    #C5
	    echo
	    echo C5_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/C5.${i}.${j}.${k} -np ~/multipool-master/C5S0/${i} ~/multipool-master/C5S7/${i}

	    #C6
	    echo
	    echo C6_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/C6.${i}.${j}.${k} -np ~/multipool-master/C6S0/${i} ~/multipool-master/C6S6/${i}

	    
	    #negative - F1
	    echo
	    echo F1_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m replicates -c ${j} -o ~/multipool-master/out/F1.${i}.${j}.${k} -np ~/multipool-master/F1/${i}
	    #negative - F2
	    echo
	    echo F2_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m replicates -c ${j} -o ~/multipool-master/out/F2.${i}.${j}.${k} -np ~/multipool-master/F2SegPool/${i}

	    #F12
	    echo
	    echo F12_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m replicates -c ${j} -o ~/multipool-master/out/F12.${i}.${j}.${k} -np ~/multipool-master/F12/${i}
	    #Combined low
	    echo
	    echo ComLow_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/ComLow.${i}.${j}.${k} -np ~/multipool-master/ComLow1/${i} ~/multipool-master/ComLow2/${i}

	    #Combined high
	    echo
	    echo ComHigh_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/ComHigh.${i}.${j}.${k} -np ~/multipool-master/ComHigh1/${i} ~/multipool-master/ComHigh2/${i}

	    #Null - contrast low
	    echo
	    echo NullLow_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/NullLow.${i}.${j}.${k} -np ~/multipool-master/NullLow1/${i} ~/multipool-master/NullLow2/${i}

	    #Null - contrast high
	    echo
	    echo NullHigh_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/NullHigh.${i}.${j}.${k} -np ~/multipool-master/NullHigh1/${i} ~/multipool-master/NullHigh2/${i}

	    #Null - contrast F2
	    echo
	    echo NullF2_${i}_${j}_${k}
	    echo
	    python ~/multipool-master/mp_inference.py -n ${k} -m contrast -c ${j} -o ~/multipool-master/out/NullF2.${i}.${j}.${k} -np ~/multipool-master/NullF21/${i} ~/multipool-master/NullF22/${i}

	done

    done

done
