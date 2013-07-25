#!/bin/bash
make
rm output_initial_forces.txt
rm config_test.cfg
echo "cutoff Epot_total Epot_local Epot_grid Epot_self force_total" >> plot.txt
for k in `seq 2 4`; do
    for l  in `seq 0 9`; do
        echo "num_iterations 1" >> config_test.cfg
        echo "num_cells 2 2 32" >> config_test.cfg
        echo "num_grids 5" >> config_test.cfg
        echo "num_gridpoints 16 16 256" >> config_test.cfg
        echo "interpolation_order 2" >> config_test.cfg
        echo "split_order 2" >> config_test.cfg
        echo "cutoff $k.$l" >> config_test.cfg
        echo "dt 5e-4" >> config_test.cfg
        config=./config_test.cfg
        ./mls_disp ./input/Ref4000atoms.input ${config} ./input/LJe0Ref4000.lammps | tee tmp.txt
        python read_forces.py >> tmp2.txt
        echo $( tail -n 1 -q tmp.txt && cat tmp2.txt) >> plot.txt
        rm tmp.txt
        rm tmp2.txt
        rm config_test.cfg
        rm output_initial_forces.txt
    done
done
