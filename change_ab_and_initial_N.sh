#!/bin/bash
START_AB=5
END_AB=100
STEP_AB=5

START_N=5
END_N=105
STEP_N=5


ab_list=()
for ((ab=START_AB; ab<=END_AB; ab=ab+STEP_AB));do
  ab_list+=($ab)
done
echo "${ab_list[@]}"

init_N_list=()
for ((N=START_N; N<=END_N; N=N+STEP_N));do
  N_list+=($N)
done
echo "${N_list[@]}"

for init_N in "${N_list[@]}";do
	echo $init_N 
	sed -i '' "s/initialN = .*/initialN = $init_N/g" variables.py 
	

	for f in "${ab_list[@]}";do
		echo $f 
		sed -i '' "s/AB_conc = .* #ug per mL/AB_conc = $f #ug per mL/g" variables.py 
		python DropTest_fct.py 
	done;
done
echo all processes complete
sed -i '' "s/AB_conc = .* #ug per mL/AB_conc = 15 #ug per mL/g" variables.py