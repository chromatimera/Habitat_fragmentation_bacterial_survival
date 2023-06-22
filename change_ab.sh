#!/bin/bash
START=15
END=75
STEP=20
python DropTest_fct.py

for ((f=START+STEP;f<=END;f=f+STEP)); do
	echo $f
	sed -i '' "s/AB_conc = .* #ug per mL/AB_conc = $f #ug per mL/g" variables.py 
	python DropTest_fct.py
done
sed -i '' "s/AB_conc = .* #ug per mL/AB_conc = 15 #ug per mL/g" variables.py