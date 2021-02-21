SHELL=/bin/bash
CORES:=$(shell nproc)
CORES43:=$(shell echo $$(($(CORES)*4/3)))
CORES_SIM:=$(shell echo $$(($(CORES)/2)))
CORES_SIM34:=$(shell echo $$(($(CORES)*3/4)))
CWD:=$(pwd)

# directory for simulations
ifndef simdir
simdir=fits
endif

ifndef area_fitscript
area_fitscript=models/joint_model.R
endif

getcores:
	echo $(CORES) $(CORES43) $(CORES_SIM) $(CORES_SIM34)

prep_sims:
	mkdir -p ${simdir}
	mkdir -p ${simdir}/log
	cp ../data/sgtf_transmission_data.rds ${simdir}
	cp ../data/targets.csv ${simdir}
	cp ../data/stp_population.csv ${simdir}
	cp ../data/i2o_rates.rds ${simdir}
	cp ${area_fitscript} ${simdir}

sims_variant_joint: prep_sims 
	Rscript gen_sim_Makefile.R --output="${simdir}/" --fitscript=$(notdir ${area_fitscript}) --variant=joint
	nice $(MAKE) -C ${simdir} -i -j$(CORES_SIM34) -l $(CORES43)

usefitstmp:
	mv ${simresults} oldfits/fits$(shell date +%Y-%m-%d.%H:%M:%S)
	mv ${simdir} ${simresults}

