#include <iostream>
#include <string>
#include <stdio.h>
#include <omp.h>

#include "mls_disp.h"
//#define UNIT_SHORTRANGE

int main(int argc, char** argv)
{
    std::string input_filename, config_filename, input_file_particles;

   if(argc<4){
      printf("Usage: %s <input_filename> <config_filename> <lammps file>\n",argv[0]);
      exit(-1);
   }else{
      input_filename = argv[1];
      config_filename = argv[2];
      input_file_particles = argv[3];
      std::cout << config_filename << std::endl;
   }

   //initialize data
   msm::Config *msm_config = new msm::Config(input_file_particles, config_filename);

#ifdef DEBUG
   msm::solver::test_verlet();
#ifdef UNIT_SHORTRANGE
   msm::ShortRange* short_range_unit = new msm::ShortRange;
   short_range_unit->computeForcesUnit();
#endif
#endif

   msm::ShortRange* short_range = new msm::ShortRange(input_filename, msm_config);
   msm::LongRange* long_range = new msm::LongRange(msm_config); 


   //start computation
   run(msm_config, short_range, long_range);

   return 0;
}
