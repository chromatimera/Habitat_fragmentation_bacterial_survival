#!/bin/bash
START=0
END=55
STEP=5

ab_list=()
for ((ab=START; ab<=END;ab=ab+STEP));do
  echo $ab
  ab_list+=($ab)
done
echo "${ab_list[@]}"


for f in "${ab_list[@]}";do
	echo $f 
	sed -i '' "s/AB_conc = .* #ug per mL/AB_conc = $f #ug per mL/g" variables.py
	python DropTest_fct.py
done
echo all processes complete
sed -i '' "s/AB_conc = .* #ug per mL/AB_conc = 15 #ug per mL/g" variables.py