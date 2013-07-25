#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <assert.h>
#include <omp.h>
#include <vector>

#include "mls_disp.h"

#define SERIES_OUTPUT

FloatType** msm::functions::gcons = 0;
FloatType** msm::functions::dgcons = 0;

#ifdef PROFILE
   double ticks_short_range = 0., ticks_long_range = 0., ticks_update = 0., direct_grid_time[] = {0., 0., 0., 0., 0.,};
#endif

/**
 * \param input_filename filename of the .lammps file which holds information about the pair coefficients
 **/
msm::Config::Config(std::string input_filename, std::string config_filename)
{
   msm::functions::init();


   using namespace std;
   //read parameters from config file
   ifstream inputfile(config_filename.c_str());
   if(inputfile) {
      string word;
      inputfile >> word;
      if(word.compare("num_iterations") == 0) inputfile >> this->num_iterations;
      inputfile >> word;
      if(word.compare("num_cells") == 0){
         inputfile >> this->num_cells[0];
         inputfile >> this->num_cells[1];
         inputfile >> this->num_cells[2];
      }
      inputfile >> word;
      if(word.compare("num_grids") == 0) inputfile >> this->num_grids;
      inputfile >> word;
      if(word.compare("num_gridpoints") == 0){
         inputfile >> this->num_gridpoints[0];
         inputfile >> this->num_gridpoints[1];
         inputfile >> this->num_gridpoints[2];
      }
      inputfile >> word;
      if(word.compare("interpolation_order") == 0) inputfile >> this->interpolation_order;
      inputfile >> word;
      if(word.compare("split_order") == 0) inputfile >> this->split_order;
      inputfile >> word;
      if(word.compare("cutoff") == 0) inputfile >> this->cutoff;
      inputfile >> word;
      if(word.compare("dt") == 0) inputfile >> this->dt;
   } else {
      cerr << "cannot open the file containing the simulation parameters " << endl;
      exit(-1);
   }

   inputfile.close();

   this->cutoff_sqr = this->cutoff * this->cutoff;
   this->stencil_size = 2 * this->interpolation_order;
   this->num_particles = 0;
   this->num_particle_types = 0;



   this->output_frequency = 1;

   vector<string> tmp_strings;
   inputfile.open(input_filename.c_str());


   FloatType boxlo[3]; 
   FloatType boxhi[3];

   //parse domain extend
   //region         box2 block 0 11.01 0 11.01 -88.08 88.08 units box 
   if(inputfile) {
      string word, word2;
      do{
         inputfile >> word;
         inputfile >> word2;
      }while(inputfile && (word.compare("region") != 0) && (word2.compare("box2") != 0));
      inputfile >> word; //skip "block"
      inputfile >> boxlo[0];
      inputfile >> boxhi[0];
      inputfile >> boxlo[1];
      inputfile >> boxhi[1];
      inputfile >> boxlo[2];
      inputfile >> boxhi[2];
#ifdef DEBUG
      printf("Domain: %f %f %f %f %f %f\n",boxlo[0],boxhi[0], boxlo[1], boxhi[1], boxlo[2], boxhi[2]);
#endif
      //make sure that the origin starts at 0,0,0
      this->domain_extend[0] = boxhi[0] - boxlo[0];
      this->domain_extend[1] = boxhi[1] - boxlo[1];
      this->domain_extend[2] = boxhi[2] - boxlo[2];
   } else {
      cerr << "cannot open file " << input_filename << endl;
      exit(-1);
   }
   this->grid_spacing[0] = this->domain_extend[0] / this->num_gridpoints[0];
   this->grid_spacing[1] = this->domain_extend[1] / this->num_gridpoints[1];
   this->grid_spacing[2] = this->domain_extend[2] / this->num_gridpoints[2];

   //reopen file
   inputfile.close();
   inputfile.open(input_filename.c_str());
   
   //read num_particles_types
   if(inputfile) {
      string word,word2;
      do{
         inputfile >> word;
      }while(inputfile && (word.compare("create_box") != 0));
      inputfile >> word2; 
      inputfile >> word; 
      if(word.compare("box2") == 0)
         this->num_particle_types = atoi(&word2[0]);
      std::cout<<"NUM TYPES: "<<num_particle_types<<std::endl;
      assert(this->num_particle_types > 0);
   } else {
      cerr << "cannot open file " << input_filename << endl;
      exit(-1);
   }

   //reopen file
   inputfile.close();
   inputfile.open(input_filename.c_str());
   
   //read num_particles
   if(inputfile) {
      string word,word2;
      do{
         inputfile >> word;
      }while(inputfile && (word.compare("create_atoms") != 0));
      inputfile >> word2; 
      inputfile >> word; 
      inputfile >> this->num_particles; 
      std::cout<<"NUM particles: "<<this->num_particles<<std::endl;
      assert(this->num_particles > 0);
   } else {
      cerr << "cannot open file " << input_filename << endl;
      exit(-1);
   }

   //allocate data
   this->epsilon = new FloatType*[this->num_particle_types];
   this->sqrt_epsilon = new FloatType[this->num_particle_types];
   this->sigma = new FloatType*[this->num_particle_types];
   this->sigma3 = new FloatType[this->num_particle_types];
   this->particle_mass = new FloatType[this->num_particle_types];
   for(int i=0;i < this->num_particle_types; ++i){
      this->epsilon[i] = new FloatType[this->num_particle_types];
      this->sigma[i] = new FloatType[this->num_particle_types];
   }
   for(int i=0;i < this->num_particle_types; ++i){
      this->particle_mass[i] = -1.0;
      for(int j=0;j < this->num_particle_types; ++j){
         this->sigma[i][j] = -1.;
         this->epsilon[i][j] = -1.;
      }
   }

   //reopen file
   inputfile.close();
   inputfile.open(input_filename.c_str());
   
   //Read sigma
   if(inputfile) {
       string word;
       do{
         inputfile >> word;
         //}while(inputfile && (word.compare("lj/cut") != 0)); // old version of input script with lj/cut
         //inputfile >> word; //skip cutoff
         //inputfile >> word;
       }while(inputfile && ((word.compare("lj/long/coul/long") != 0) && (word.compare("lj/cut") != 0)));
      if(word.compare("lj/long/coul/long") == 0){
          inputfile >> word; // long
          inputfile >> word; // off
          inputfile >> word; // value of cutoff
      }
      else if(word.compare("lj/cut") == 0){
          inputfile >> word; //skip cutoff
      }
      inputfile >> word;
      while(word.compare("pair_coeff") == 0){
         int i,j;
         inputfile >> i;
         inputfile >> j;
         --i; --j;
         inputfile >> this->epsilon[i][j];
         inputfile >> this->sigma[i][j];
         this->epsilon[j][i] = this->epsilon[i][j];
         this->sigma[j][i] = this->sigma[i][j];
         inputfile >> word;
      }
   } else {
      cerr << "cannot open file " << input_filename << endl;
      exit(-1);
   }

   //reopen file
   inputfile.close();
   inputfile.open(input_filename.c_str());
   
   //Read epsilon 
   if(inputfile) {
      string word;
      do{
         inputfile >> word;
      }while(inputfile && (word.compare("mass") != 0));

      while(word.compare("mass") == 0){
         int i;
         inputfile >> i;
         --i; 
         inputfile >> this->particle_mass[i];
         inputfile >> word;
      }
   } else {
      cerr << "cannot open file " << input_filename << endl;
      exit(-1);
   }
   
   //test sigma and epsilon
   for(int i=0;i < this->num_particle_types; ++i){
      if(this->particle_mass[i] <= 0.0){
          std::cout<< "ERROR: particle mass <= 0.0!\n";
            exit(-1);
      }
      for(int j=0;j < this->num_particle_types; ++j){
#ifdef DEBUG
         printf("e: %f s: %f\n",this->epsilon[j][i],this->sigma[j][i]);
#endif
         if(this->epsilon[i][j] < 0 || this->sigma[i][j] < 0){
            std::cout<< "ERROR: not all sigma and epsilon read!\n";
            exit(-1);
         }
      }
   }
   for(int i=0;i < this->num_particle_types; ++i){
       this->sigma3[i]= sqr(this->sigma[i][i]) * this->sigma[i][i];
       this->sqrt_epsilon[i]= sqrt(this->epsilon[i][i]);
   }

   assert(this->cutoff > 0 && this->cutoff_sqr > 0 && this->num_grids > 0);
   assert(this->split_order <= 6 && this->split_order >= 2);
   assert(this->interpolation_order <= 5 && this->interpolation_order >= 2);
}

msm::Config::~Config()
{
   //fee data
   for(int i=0;i < this->num_particle_types; ++i){
      delete[] this->epsilon[i];
      delete[] this->sigma[i];
   }
   delete[] this->particle_mass;
   delete[] this->epsilon;
   delete[] this->sigma;
}

/*-----------------------------------------------------------------------------
 *  namespace:solver 
 *-----------------------------------------------------------------------------*/
    /**
     * Integrate All given particles in time using the velocity verlet scheme.
     *
     * \param[in]     dt                time step    
     * \param[in]     num_particles     number of particles to integrate
     * \param[in]     particleForce     this array stores the force acting on each particle \
     *                                  at the current time step.
     * \param[in]     particleForceOld  this array stores the force acting on each particle \
     *                                  at the previous time step.
     * \param[in,out] particleVel       this array stores the velocity of each particle.
     * \param[in,out] particlePos       this array stores the position of each particle.
     **/
    void msm::solver::verlet( const FloatType dt,
            const int num_particles,
            int* particle_type,
            FloatType* particle_mass,
            FloatType** particle_force,
            FloatType** particle_force_old, 
            FloatType** particle_vel, 
            FloatType** particle_pos)
    {
       int i,d;
       for(i=0;i < num_particles ;++i){

          FloatType accel = dt * .5 / particle_mass[particle_type[i]]; 


          for(d = 0 ; d < 3; ++d){
             FloatType update = dt * (particle_vel[i][d] 
                   + accel * particle_force[i][d]);

             particle_pos[i][d] += update;
             particle_vel[i][d] += accel * (particle_force_old[i][d] 
                   + particle_force[i][d]);
          }
       }
    }

void msm::solver::test_verlet()
{
   int num_particles = 1024;
   int num_particle_types = 10;
   FloatType dt = 0.0005;

   int* particle_type;
   FloatType* particle_mass;
   FloatType** particle_force;
   FloatType** particle_force_old;
   FloatType** particle_vel; 
   FloatType** particle_vel_tmp; 
   FloatType** particle_pos;
   FloatType** particle_pos_tmp;

   memory::create(particle_pos, num_particles, 3);
   memory::create(particle_pos_tmp, num_particles, 3);
   memory::create(particle_vel, num_particles, 3);
   memory::create(particle_vel_tmp, num_particles, 3);
   memory::create(particle_force, num_particles, 3);
   memory::create(particle_force_old, num_particles, 3);
   particle_mass             = new FloatType[num_particle_types];
   particle_type             = new int[num_particles];

   //test for each direction
   for(int d=0;d < 3; ++d){

      for(int i=0;i < num_particle_types;++i){
         particle_mass[i] = rand()/((FloatType)RAND_MAX);
      }

      for(int i=0;i < num_particles ;++i){
         particle_type[i] = rand()%num_particle_types;
         particle_pos[i][0] = rand()/((FloatType)RAND_MAX);
         particle_pos[i][1] = rand()/((FloatType)RAND_MAX);
         particle_pos[i][2] = rand()/((FloatType)RAND_MAX);
         particle_pos_tmp[i][0] = particle_pos[i][0];
         particle_pos_tmp[i][1] = particle_pos[i][1];
         particle_pos_tmp[i][2] = particle_pos[i][2];
         particle_vel[i][0] = 0.0;
         particle_vel[i][1] = 0.0;
         particle_vel[i][2] = 0.0;
         particle_vel_tmp[i][0] = particle_vel[i][0];
         particle_vel_tmp[i][1] = particle_vel[i][1];
         particle_vel_tmp[i][2] = particle_vel[i][2];
      }

      //force in x dim must not affect y and z positions/velocity, some goes for y and z
      for(int i=0;i < num_particles ;++i){
         particle_force[i][0] = (d != 0) ? 0.0 : rand()/((FloatType)RAND_MAX);
         particle_force[i][1] = (d != 1) ? 0.0 : rand()/((FloatType)RAND_MAX);
         particle_force[i][2] = (d != 2) ? 0.0 : rand()/((FloatType)RAND_MAX);
         particle_force_old[i][0] = (d != 0) ? 0.0 : rand()/((FloatType)RAND_MAX);
         particle_force_old[i][1] = (d != 1) ? 0.0 : rand()/((FloatType)RAND_MAX);
         particle_force_old[i][2] = (d != 2) ? 0.0 : rand()/((FloatType)RAND_MAX);
      }

      //three steps of verlet method in order to change the velocity over time
      msm::solver::verlet( dt, num_particles, particle_type,
            particle_mass, particle_force, particle_force_old, 
            particle_vel, particle_pos);
      msm::solver::verlet( dt, num_particles, particle_type,
            particle_mass, particle_force, particle_force_old, 
            particle_vel, particle_pos);
      msm::solver::verlet( dt, num_particles, particle_type,
            particle_mass, particle_force, particle_force_old, 
            particle_vel, particle_pos);

      //force in x dim must not affect y and z positions/velocity
      for(int i=0;i < num_particles ;++i){
         if(particle_pos_tmp[i][(d+1)%3] - particle_pos[i][(d+1)%3] != 0.0 ||
               particle_pos_tmp[i][(d+2)%3] - particle_pos[i][(d+2)%3] != 0.0 ||
               particle_vel_tmp[i][(d+1)%3] - particle_vel[i][(d+1)%3] != 0.0 ||
               particle_vel_tmp[i][(d+2)%3] - particle_vel[i][(d+2)%3] != 0.0){
            std::cerr<<"ERROR: verlet() seems to be buggy!\n";
            exit(-1);
         }
      }
   }

   free(&(particle_vel[0][0]));
   free(&(particle_vel_tmp[0][0]));
   free(&(particle_pos[0][0]));
   free(&(particle_pos_tmp[0][0]));
   free(&(particle_force[0][0]));
   free(&(particle_force_old[0][0]));
   free(particle_mass);
   free(particle_type);

   std::cout<<"VERLET TEST SUCCESSFULL!\n";
}

/*-----------------------------------------------------------------------------
 *  Class: LongRange
 *-----------------------------------------------------------------------------*/
msm::LongRange::LongRange(msm::Config* config)
    {
        this->num_grids = config->num_grids;
        if(num_grids<2){
            std::cerr<<"ERROR: number of grids must be greater or equal to 2. If you require less grids please contact us.\n";
            exit(-1);
        }
        this->potential_energy = 0.0;
        int num_gridpoints[3];
        //copy values in order to leave config untouched
        num_gridpoints[0] = config->num_gridpoints[0];
        num_gridpoints[1] = config->num_gridpoints[1];
        num_gridpoints[2] = config->num_gridpoints[2];
        this->config = config;
        this->grids = new Grid*[this->num_grids];
        this->grids[0] = new FinestGrid;
        this->grids[0]->init(config, 0);

        for(int i=1;i < num_grids-1 ;++i){
            num_gridpoints[0] /= 2; num_gridpoints[1] /= 2; num_gridpoints[2] /= 2;
            if(!(num_gridpoints[0] > 0 && num_gridpoints[1] > 0 && num_gridpoints[2] > 0)){
                std::cout << "The number of points on a grid became zero. Adapt the number of points on the finest grid or the number of grids" << std::endl;
                exit(-1);
            }
            this->grids[i] = new Grid;
            this->grids[i]->init(config, i);
        }

        num_gridpoints[0] /= 2; num_gridpoints[1] /= 2; num_gridpoints[2] /= 2;
        if (!((num_gridpoints[0] > 0) && (num_gridpoints[1] > 0) && (num_gridpoints[2] > 0))){
                std::cout << "The number of points on a grid became zero. Adapt the number of points on the finest grid or the number of grids" << std::endl;
                exit(-1);
            }
        this->grids[this->num_grids-1] = new CoarsestGrid;
        this->grids[this->num_grids-1]->init(config, num_grids-1);

        //precompute phi for restriction and prolongation
        phi_grid= new FloatType[2 * config->interpolation_order+1];
        phi_grid_index = new int[2 * config->interpolation_order+1];
        if(config->interpolation_order == 2){
           phi_grid[0] = functions::phi(-1.5, 2 * config->interpolation_order);
           phi_grid_index[0] = -3;
           phi_grid[1] = functions::phi(-0.5, 2 * config->interpolation_order);
           phi_grid_index[1] = -1;
           phi_grid[2] = functions::phi(0., 2 * config->interpolation_order);
           phi_grid_index[2] = 0;
           phi_grid[3] = functions::phi(0.5, 2 * config->interpolation_order);
           phi_grid_index[3] = 1;
           phi_grid[4] = functions::phi(1.5, 2 * config->interpolation_order);
           phi_grid_index[4] = 3;
        }else if (config->interpolation_order == 3){
           phi_grid[0] = functions::phi(-2.5, 2 * config->interpolation_order);
           phi_grid_index[0] = -5;
           phi_grid[1] = functions::phi(-1.5, 2 * config->interpolation_order);
           phi_grid_index[1] = -3;
           phi_grid[2] = functions::phi(-0.5, 2 * config->interpolation_order);
           phi_grid_index[2] = -1;
           phi_grid[3] = functions::phi(0., 2 * config->interpolation_order);
           phi_grid_index[3] = 0;
           phi_grid[4] = functions::phi(0.5, 2 * config->interpolation_order);
           phi_grid_index[4] = 1;
           phi_grid[5] = functions::phi(1.5, 2 * config->interpolation_order);
           phi_grid_index[5] = 3;
           phi_grid[6] = functions::phi(2.5, 2 * config->interpolation_order);
           phi_grid_index[6] = 5;
        }else if (config->interpolation_order == 4){
           phi_grid[0] = functions::phi(-3.5, 2 * config->interpolation_order);
           phi_grid_index[0] = -7;
           phi_grid[1] = functions::phi(-2.5, 2 * config->interpolation_order);
           phi_grid_index[1] = -5;
           phi_grid[2] = functions::phi(-1.5, 2 * config->interpolation_order);
           phi_grid_index[2] = -3;
           phi_grid[3] = functions::phi(-0.5, 2 * config->interpolation_order);
           phi_grid_index[3] = -1;
           phi_grid[4] = functions::phi(0., 2 * config->interpolation_order);
           phi_grid_index[4] = 0;
           phi_grid[5] = functions::phi(0.5, 2 * config->interpolation_order);
           phi_grid_index[5] = 1;
           phi_grid[6] = functions::phi(1.5, 2 * config->interpolation_order);
           phi_grid_index[6] = 3;
           phi_grid[7] = functions::phi(2.5, 2 * config->interpolation_order);
           phi_grid_index[7] = 5;
           phi_grid[8] = functions::phi(3.5, 2 * config->interpolation_order);
           phi_grid_index[8] = 7;
        }else if (config->interpolation_order == 5){
           phi_grid[0] = functions::phi(-4.5, 2 * config->interpolation_order);
           phi_grid_index[0] = -9;
           phi_grid[1] = functions::phi(-3.5, 2 * config->interpolation_order);
           phi_grid_index[1] = -7;
           phi_grid[2] = functions::phi(-2.5, 2 * config->interpolation_order);
           phi_grid_index[2] = -5;
           phi_grid[3] = functions::phi(-1.5, 2 * config->interpolation_order);
           phi_grid_index[3] = -3;
           phi_grid[4] = functions::phi(-0.5, 2 * config->interpolation_order);
           phi_grid_index[4] = -1;
           phi_grid[5] = functions::phi(0., 2 * config->interpolation_order);
           phi_grid_index[5] = 0;
           phi_grid[6] = functions::phi(0.5, 2 * config->interpolation_order);
           phi_grid_index[6] = 1;
           phi_grid[7] = functions::phi(1.5, 2 * config->interpolation_order);
           phi_grid_index[7] = 3;
           phi_grid[8] = functions::phi(2.5, 2 * config->interpolation_order);
           phi_grid_index[8] = 5;
           phi_grid[9] = functions::phi(3.5, 2 * config->interpolation_order);
           phi_grid_index[9] = 7;
           phi_grid[10] = functions::phi(4.5, 2 * config->interpolation_order);
           phi_grid_index[10] = 9;
        }else{
           std::cerr<< "Error: Interpolation order not yet supported.\n";
        }

        {
           //precompute gamma_grid for direct() functions, except for last grid
           //this gamma grid is computed for the finest grid. However, other grids can use
           //the same values and multiply them with a proper factor
           FloatType cutoff2_sqr = sqr(2 *config->cutoff);
           //For finest grid, the cutoff is equal to 2*config->cutoff.
           //number of gridpoints within the stencil in each direction
           int num_points_cut_x = (2 * config->cutoff)/config->grid_spacing[0];
           int num_points_cut_y = (2 * config->cutoff)/config->grid_spacing[1];
           int num_points_cut_z = (2 * config->cutoff)/config->grid_spacing[2];
           FloatType grid_spacing_x_sqr = sqr(config->grid_spacing[0]);
           FloatType grid_spacing_y_sqr = sqr(config->grid_spacing[1]);
           FloatType grid_spacing_z_sqr = sqr(config->grid_spacing[2]);
           //count number of gridpoints inside of the stencil (depending on cutoff)
           int counter = 0;
           for(int i= -num_points_cut_z;i <= num_points_cut_z; ++i){
              for(int j= -num_points_cut_y;j <= num_points_cut_y; ++j){
                 for(int k= -num_points_cut_x;k <= num_points_cut_x; ++k){
                    FloatType distance_sqr = sqr(i) * grid_spacing_z_sqr + sqr(j) * grid_spacing_y_sqr + sqr(k) * grid_spacing_x_sqr;
                    if(distance_sqr < cutoff2_sqr)
                       ++counter;
                 }
              }
           }
           assert(counter <= (2*num_points_cut_x+1) * (2*num_points_cut_y+1) * (2*num_points_cut_z+1));
           printf("Number of gridpoints within cutoff for the finest grid: %d\n", counter);

           gamma_grid = new FloatType[counter];
           gamma_grid_index[0] = new int[counter];
           gamma_grid_index[1] = new int[counter];
           gamma_grid_index[2] = new int[counter];
           gamma_grid_size = counter;

           counter = 0;
           FloatType cutoff_sqr = sqr(config->cutoff);
           for(int i= -num_points_cut_z;i <= num_points_cut_z; ++i){
              for(int j= -num_points_cut_y;j <= num_points_cut_y; ++j){
                 for(int k= -num_points_cut_x;k <= num_points_cut_x; ++k){
                    FloatType distance_sqr = sqr(i) * grid_spacing_z_sqr + sqr(j) * grid_spacing_y_sqr + sqr(k) * grid_spacing_x_sqr;
                    if(distance_sqr < cutoff2_sqr){
                       FloatType tmp = 0.0;
                       if(distance_sqr <= cutoff_sqr){
                           tmp = msm::functions::gamma(distance_sqr/cutoff_sqr, config->split_order);
                       }else{
                           tmp = distance_sqr/cutoff_sqr; // (r/cutoff)^2
                           tmp = 1.0 / (sqr(tmp) * tmp); // (r/cutoff)^-6
                       }
                       //1/cutoff^6 * tmp - (2*cutoff)^-6 * gamma((r/(2*cutoff))^2)
                       gamma_grid[counter] = 1.0 / (sqr(cutoff_sqr) * cutoff_sqr) * tmp 
                          - 1.0 / (sqr(cutoff2_sqr) * cutoff2_sqr) * msm::functions::gamma(distance_sqr/cutoff2_sqr, config->split_order);
                       gamma_grid_index[0][counter] = k;
                       gamma_grid_index[1][counter] = j;
                       gamma_grid_index[2][counter] = i;
                       ++counter;
                    }
                 }
              }
           }
           assert(counter == gamma_grid_size && counter > 0);
        }
    }

    msm::LongRange::~LongRange()
    {
        for(int i=0;i < this->num_grids ;++i){
            delete(grids[i]);
        }
        delete [] this->grids;
        delete [] phi_grid;
        delete [] phi_grid_index;

        delete [] gamma_grid;
        delete [] gamma_grid_index[0];
        delete [] gamma_grid_index[1];
        delete [] gamma_grid_index[2];
    }

    /**
     * This function computes one step of the LongRange contribution to the final force on multiple grids.
     * After this step all the forces of all particles in every cell will be updated.
     * 
     * \param[in]       num_cells   number of cells in each dimension. num_cells[0] : cells in x-dim, ... 
     * \param[in,out]   cells       3D pointer to Cells. Cell[z-dim][y-dim][x-dim].
     **/
    void msm::LongRange::computeForces(int num_cells[3], msm::Cell*** cells)
    {

        this->potential_energy = 0.0;
        //map particle xi to grid
        msm::FinestGrid* finest = dynamic_cast<msm::FinestGrid*>(this->grids[0]);
        assert(finest);

        finest->resetGrid();
        finest->anterpolate(cells, num_cells, this->config);

        //propagate xi to lower level grids
        for(int i=0;i < num_grids-1 ;++i){
           grids[i+1]->resetGrid(); 
           grids[i]->restriction(this->grids[i+1], this->config->stencil_size, this->phi_grid, this->phi_grid_index);
#ifdef PROFILE
        double start_time = omp_get_wtime();
#endif
           grids[i]->direct(this->config->cutoff, this->gamma_grid, this->gamma_grid_index, this->gamma_grid_size);
#ifdef PROFILE
        direct_grid_time[i] += omp_get_wtime()-start_time;
#endif
        }
        //coarsest grid
#ifdef PROFILE
        double start_time = omp_get_wtime();
#endif
            grids[num_grids-1]->direct(this->config->cutoff);
#ifdef PROFILE
        direct_grid_time[num_grids-1] += omp_get_wtime()-start_time;
#endif

        //propagate higher level contribution to lower level grids
        for(int i=num_grids-1 ;i > 0 ;--i){
            grids[i]->prolongation(this->grids[i-1], this->config->stencil_size, this->phi_grid, this->phi_grid_index);
        }

        //map grid xi to particles
        this->potential_energy += finest->interpolate(cells, num_cells, this->config);
    }

    /*-----------------------------------------------------------------------------
     *  Class: Grid
     *-----------------------------------------------------------------------------*/
    msm::Grid::Grid(){
        this->spacing_grid[0]          = -1.0;
        this->spacing_grid[1]          = -1.0;
        this->spacing_grid[2]          = -1.0;
        this->interpolation_order = -1;
        this->level               = -1;
        this->gridpoint_potential = NULL;
        this->gridpoint_force_x = NULL;
        this->gridpoint_force_y = NULL;
        this->gridpoint_force_z = NULL;
        this->gridpoint_xi    = NULL;
        this->direct_stencil   = NULL;
    }

    msm::Grid::~Grid(){
       free3DArray(this->gridpoint_potential, num_gridpoints[1], num_gridpoints[2]);
       if(this->gridpoint_force_x)
           free3DArray(this->gridpoint_force_x, num_gridpoints[1], num_gridpoints[2]);
       if(this->gridpoint_force_y)
           free3DArray(this->gridpoint_force_y, num_gridpoints[1], num_gridpoints[2]);
       if(this->gridpoint_force_z)
           free3DArray(this->gridpoint_force_z, num_gridpoints[1], num_gridpoints[2]);
       free3DArray(this->gridpoint_xi, num_gridpoints[1], num_gridpoints[2]);

       if(direct_stencil)
         free3DArray(this->direct_stencil, num_gridpoints[1], num_gridpoints[2]);
    }

    void msm::Grid::init(msm::Config* config, int level)
    {
        //Grid spacing is only used by finest and coarsest grid (see coarsestGrid::init)
        this->spacing_grid[0] = config->grid_spacing[0];
        this->spacing_grid[1] = config->grid_spacing[1];
        this->spacing_grid[2] = config->grid_spacing[2];

        assert(config->num_gridpoints[0] % (1<<level) == 0);
        assert(config->num_gridpoints[1] % (1<<level) == 0);
        assert(config->num_gridpoints[2] % (1<<level) == 0);
        this->num_gridpoints[0] = config->num_gridpoints[0]/(1<<level);
        this->num_gridpoints[1] = config->num_gridpoints[1]/(1<<level);
        this->num_gridpoints[2] = config->num_gridpoints[2]/(1<<level);

        this->interpolation_order = config->interpolation_order;
        this->level               = level;

        allocate3DArray<FloatType> (this->gridpoint_potential, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        allocate3DArray<FloatType> (this->gridpoint_xi, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
       
        this->direct_aux = 1.0/(pow(64.0, this->level)); 
        
        // cutoff on grid is 2 << level * a (grid 0: cutoff 2 )
        // cutoff on the grid divided by domain_extend + 1 = how often periodic boundaries can occur = how often one has to add num_gridpoints to get a positive number 
        this->direct_periodic_int_mult_of_num_gridpoints[0] = static_cast<int>(floor( ((2<<level)*config->cutoff) / config->domain_extend[0]) + 1) * this->num_gridpoints[0];
        this->direct_periodic_int_mult_of_num_gridpoints[1] = static_cast<int>(floor( ((2<<level)*config->cutoff) / config->domain_extend[1]) + 1) * this->num_gridpoints[1];
        this->direct_periodic_int_mult_of_num_gridpoints[2] = static_cast<int>(floor( ((2<<level)*config->cutoff) / config->domain_extend[2]) + 1) * this->num_gridpoints[2];
    }

    void msm::Grid::direct(FloatType cutoff, FloatType *gamma_grid, int **gamma_grid_index, int gamma_grid_size){

        assert(this->gridpoint_potential[0][0][0] == 0);

        //iterate over all gridpoints
        for(int l = 0; l < gamma_grid_size; ++l){
           FloatType tmp = this->direct_aux * gamma_grid[l];
           for(int i=0;i < this->num_gridpoints[2] ;++i){
              int ii = (gamma_grid_index[2][l] + i + this->direct_periodic_int_mult_of_num_gridpoints[2])%this->num_gridpoints[2]; 
              for(int j=0;j < this->num_gridpoints[1] ;++j){
                 int jj = (gamma_grid_index[1][l] + j + this->direct_periodic_int_mult_of_num_gridpoints[1])%this->num_gridpoints[1];
                 for(int k=0;k < this->num_gridpoints[0] ;++k){
                    //periodic boundaries
                    int kk = (gamma_grid_index[0][l] + k + this->direct_periodic_int_mult_of_num_gridpoints[0])%this->num_gridpoints[0];
                    this->gridpoint_potential[i][j][k] += tmp * this->gridpoint_xi[ii][jj][kk];
                 }
              }
           }
        }
    }

    void msm::Grid::resetGrid()
    {
        for(int i=0;i < this->num_gridpoints[2] ;++i){
            for(int j=0;j < this->num_gridpoints[1] ;++j){
                for(int k=0;k < this->num_gridpoints[0] ;++k){
                    this->gridpoint_xi[i][j][k] = 0.0;
                    this->gridpoint_potential[i][j][k] = 0.0;
                }
            }
        }
    }

    /**
     * \param phi_grid array of size 2 * interpolation order + 1. Phi-function evaluated at neighboring gridpoints (zeros are not stored).
     * \param phi_grid_inex corresponding indecies to phi_grid
     */
    void msm::Grid::restriction(Grid* dst_grid, const int stencil_size, FloatType *phi_grid, int* phi_grid_index)
    {
        FloatType*** dst_xi = dst_grid->getGridpointXi();

        //loop over all gridpoints of the higher-level (coarser) grid and update their charges
        for(int i=0; i < dst_grid->getNumGridpoints()[2] ;++i){
            for(int j=0; j < dst_grid->getNumGridpoints()[1] ;++j){
                for(int k=0; k < dst_grid->getNumGridpoints()[0] ;++k){

                    assert(dst_xi[i][j][k] == 0.0);
                    //determine indecies of gridpoint belonging to the finer grid which alignes
                    //with the current index of the coarser grid
                    int fi = 2*i, fj = 2*j, fk = 2*k;

                    //WARNING: 2*interpolation_order must be equal to stencil_size
                    //loop over stencil of the finer grid and update the charge of the coarser grid
                    for(int ii = 0; ii < stencil_size + 1; ++ii){
                       // multiplication with stencil_size is necesarry because of interpolation order and coarseness of the grid
                       int tmp_ii = (fi + phi_grid_index[ii] + stencil_size * this->num_gridpoints[2]) % this->num_gridpoints[2];
                       for(int jj = 0; jj < stencil_size + 1; ++jj){
                           int tmp_jj = (fj + phi_grid_index[jj] + stencil_size * this->num_gridpoints[1]) % this->num_gridpoints[1];
                           for(int kk = 0; kk < stencil_size + 1; ++kk){
                              int tmp_kk = (fk + phi_grid_index[kk] + stencil_size * this->num_gridpoints[0]) % this->num_gridpoints[0];

                              assert(ii >= 0 && ii < stencil_size + 1 &&
                                    jj >= 0 && jj < stencil_size + 1 &&
                                    kk >= 0 && kk < stencil_size + 1);
                              assert( tmp_ii >= 0 && tmp_ii < this->num_gridpoints[2] && //out of bounds check
                                    tmp_jj >= 0 && tmp_jj < this->num_gridpoints[1] &&
                                    tmp_kk >= 0 && tmp_kk < this->num_gridpoints[0]); 

                              dst_xi[i][j][k] += phi_grid[ii] * phi_grid[jj] * phi_grid[kk] 
                                 * this->gridpoint_xi[tmp_ii][tmp_jj][tmp_kk];
                           }
                       }
                    }

                }
            }
        }
    }

void msm::Grid::prolongation(Grid* dst_grid, const int stencil_size, FloatType* phi_grid, int* phi_grid_index)
    {
        FloatType*** dst_potential = dst_grid->getGridpointPotential();
        
        int* dst_num_gridpoints = dst_grid->getNumGridpoints();

        /* 
         * loop over all gridpoints of this grid (coarser) and compute the 
         * contribution of the charges for the lower-level grid 
         * */
        for(int i=0;i < this->num_gridpoints[2] ;++i){
            for(int j=0;j < this->num_gridpoints[1] ;++j){
                for(int k=0;k < this->num_gridpoints[0] ;++k){

                    //determine central index of the lower-level grid
                    int fi = 2*i, fj = 2*j, fk = 2*k;

                    //WARNING: 2*interpolation_order must be equal to stencil_size
                    //loop over the stencil and update the potential of the destination grid
                    for(int ii = 0; ii < stencil_size + 1; ++ii){
                       int tmp_ii = (fi + phi_grid_index[ii] + stencil_size * dst_num_gridpoints[2]) % dst_num_gridpoints[2];
                       for(int jj = 0; jj < stencil_size + 1; ++jj){
                           int tmp_jj = (fj + phi_grid_index[jj] + stencil_size * dst_num_gridpoints[1]) % dst_num_gridpoints[1];
                           for(int kk = 0; kk < stencil_size + 1; ++kk){
                              int tmp_kk = (fk + phi_grid_index[kk] + stencil_size * dst_num_gridpoints[0]) % dst_num_gridpoints[0];

                              assert(ii >= 0 && ii < stencil_size + 1 &&
                                    jj >= 0 && jj < stencil_size + 1 &&
                                    kk >= 0 && kk < stencil_size + 1);
                              assert( tmp_ii >= 0 && tmp_ii < dst_num_gridpoints[2] && //out of bounds check
                                    tmp_jj >= 0 && tmp_jj < dst_num_gridpoints[1] &&
                                    tmp_kk >= 0 && tmp_kk < dst_num_gridpoints[0]); 

                               dst_potential[tmp_ii][tmp_jj][tmp_kk] +=  phi_grid[ii] * phi_grid[jj] * phi_grid[kk] 
                                  * this->gridpoint_potential[i][j][k];
                           }
                       }
                    }
                }
            }
        }
    }
                       
    /**
     * \param[in]  particle_pos  pointer to all particles in a current cell. \
                                 Dev: this is a little overkill but avoids unnecessary memory copies.
     * \param[in]  particle_idx  particle index within particle_pos.
     * \param[out] lower_grid_idx stores the lower, left, front gridpoint index that is affected by the given particle*/
    void msm::FinestGrid::find_lowest_gridpoint_idx(FloatType** particle_pos, const int particle_idx, int (&lower_grid_idx)[3])
    {
       lower_grid_idx[0] = static_cast<int>(floor(particle_pos[particle_idx][0] / this->spacing_grid[0])) - (this->interpolation_order - 1); //x-dim
       lower_grid_idx[1] = static_cast<int>(floor(particle_pos[particle_idx][1] / this->spacing_grid[1])) - (this->interpolation_order - 1); //y-dim
       lower_grid_idx[2] = static_cast<int>(floor(particle_pos[particle_idx][2] / this->spacing_grid[2])) - (this->interpolation_order - 1); //z-dim
    }

    /**
     * \param[in] cells 3D    pointer to all cells. 
     * \param[in] num_cells   number of cells in each dimension. num_cells[0]=x-dim;num_cells[1]=y-dim;...
     **/
    void msm::FinestGrid::anterpolate(msm::Cell*** cells, const int num_cells[3], msm::Config* config)
    {
        assert(this->gridpoint_xi[0][0][0] == 0);

        FloatType* phi[3];
        for(int ii = 0; ii < 3; ++ii){
           phi[ii] = new FloatType[config->stencil_size];
        }

        //iterate over all cells
        for(int i=0;i < num_cells[2] ;++i){
            for(int j=0;j < num_cells[1] ;++j){
                for(int k=0;k < num_cells[0] ;++k){

                    msm::Cell&        current_cell   = cells[i][j][k];
                    FloatType**   particle_pos   = current_cell.getParticlePos();
                    int*         particle_type  = current_cell.getParticleType();

                    //iterate over all particles of a cell
                    int num_particles = current_cell.getNumParticles();
                    for(int particle_idx=0; particle_idx < num_particles; ++particle_idx){
                       int lower_grid_idx[3];
                       find_lowest_gridpoint_idx(particle_pos, particle_idx, lower_grid_idx);
                       
                       //compute phi
                       for(int ii = 0; ii < config->stencil_size; ++ii){
                           phi[0][ii] = msm::functions::phi(lower_grid_idx[0] - (particle_pos[particle_idx][0] / this->spacing_grid[0]) + ii, 2 * config->interpolation_order);
                           phi[1][ii] = msm::functions::phi(lower_grid_idx[1] - (particle_pos[particle_idx][1] / this->spacing_grid[1]) + ii, 2 * config->interpolation_order);
                           phi[2][ii] = msm::functions::phi(lower_grid_idx[2] - (particle_pos[particle_idx][2] / this->spacing_grid[2]) + ii, 2 * config->interpolation_order);
                       }

                       //loop over all grid points (small stencil) which are affected by particle_idx
                       for(int ii = 0; ii < config->stencil_size; ++ii){
                           for(int jj = 0; jj < config->stencil_size; ++jj){
                               for(int kk = 0; kk < config->stencil_size; ++kk){

                                   //periodic boundary
                                   int idx_z = (lower_grid_idx[2] + ii + num_gridpoints[2]) % num_gridpoints[2];
                                   int idx_y = (lower_grid_idx[1] + jj + num_gridpoints[1]) % num_gridpoints[1];
                                   int idx_x = (lower_grid_idx[0] + kk + num_gridpoints[0]) % num_gridpoints[0];

                                   assert( idx_z >= 0 && idx_z < this->num_gridpoints[2] && //out of bounds check
                                           idx_y >= 0 && idx_y < this->num_gridpoints[1] &&
                                           idx_x >= 0 && idx_x < this->num_gridpoints[0]); 

                                   this->gridpoint_xi[idx_z][idx_y][idx_x] += 
                                       phi[0][kk] * phi[1][jj] * phi[2][ii] * 
                                       config->sigma3[particle_type[particle_idx]] * 
                                       config->sqrt_epsilon[particle_type[particle_idx]];
                               }
                           }
                       } 
                    }//end particle_idx
                }
            }
        }

        for(int ii = 0; ii < 3; ++ii){
           delete[] phi[ii];
        }
    }

   /**
    * \return Potential from the long-range part
    */
FloatType msm::FinestGrid::interpolate(msm::Cell*** cells, const int num_cells[3], msm::Config* config)
    {
        FloatType* phi[3];
        FloatType* dphi[3];
        for(int ii = 0; ii < 3; ++ii){
           phi[ii] = new FloatType[config->stencil_size];
           dphi[ii] = new FloatType[config->stencil_size];
        }

        FloatType potential = 0.0;

        //iterate over all cells
        for(int i=0;i < num_cells[2] ;++i){
            for(int j=0;j < num_cells[1] ;++j){
                for(int k=0;k < num_cells[0] ;++k){

                    msm::Cell&    current_cell  = cells[i][j][k];
                    FloatType**      particle_pos  = current_cell.getParticlePos();
                    FloatType**      particle_force_grid = current_cell.getParticleForceGrid();
                    int*            particle_type = current_cell.getParticleType();

                    //iterate over all particles of a cell
                    int num_particles = current_cell.getNumParticles();
                    for(int particle_idx=0; particle_idx < num_particles; ++particle_idx){
                       int lower_grid_idx[3];
                       find_lowest_gridpoint_idx(particle_pos, particle_idx, lower_grid_idx);
                       
                       //compute phi and its derivative dphi
                       for(int ii = 0; ii < config->stencil_size; ++ii){
                           phi[0][ii] = msm::functions::phi(lower_grid_idx[0] - (particle_pos[particle_idx][0] / this->spacing_grid[0]) + ii, 2 * config->interpolation_order);
                           phi[1][ii] = msm::functions::phi(lower_grid_idx[1] - (particle_pos[particle_idx][1] / this->spacing_grid[1]) + ii, 2 * config->interpolation_order);
                           phi[2][ii] = msm::functions::phi(lower_grid_idx[2] - (particle_pos[particle_idx][2] / this->spacing_grid[2]) + ii, 2 * config->interpolation_order);

                           dphi[0][ii] = msm::functions::dphi(lower_grid_idx[0] - (particle_pos[particle_idx][0] / this->spacing_grid[0]) + ii, 2 * config->interpolation_order);
                           dphi[1][ii] = msm::functions::dphi(lower_grid_idx[1] - (particle_pos[particle_idx][1] / this->spacing_grid[1]) + ii, 2 * config->interpolation_order);
                           dphi[2][ii] = msm::functions::dphi(lower_grid_idx[2] - (particle_pos[particle_idx][2] / this->spacing_grid[2]) + ii, 2 * config->interpolation_order);
                       }

                       FloatType force[3] = {0.0, 0.0, 0.0};
                       //loop over all grid points (small stencil) which contribute to particle_idx
                       FloatType tmp_potential = 0.0;
                       for(int ii = 0; ii < config->stencil_size; ++ii){
                           for(int jj = 0; jj < config->stencil_size; ++jj){
                               for(int kk = 0; kk < config->stencil_size; ++kk){

                                  //periodic boundary
                                  int idx_z = (lower_grid_idx[2] + ii + num_gridpoints[2]) % num_gridpoints[2];
                                  int idx_y = (lower_grid_idx[1] + jj + num_gridpoints[1]) % num_gridpoints[1];
                                  int idx_x = (lower_grid_idx[0] + kk + num_gridpoints[0]) % num_gridpoints[0];

                                  assert( idx_z >= 0 && idx_z < this->num_gridpoints[2] && //out of bounds check
                                        idx_y >= 0 && idx_y < this->num_gridpoints[1] &&
                                        idx_x >= 0 && idx_x < this->num_gridpoints[0]); 

                                  tmp_potential += phi[0][kk] * phi[1][jj] * phi[2][ii] * gridpoint_potential[idx_z][idx_y][idx_x];
                                  force[0] += dphi[0][kk] / this->spacing_grid[0] *  phi[1][jj] *  phi[2][ii] * gridpoint_potential[idx_z][idx_y][idx_x]; 
                                  force[1] +=  phi[0][kk] * dphi[1][jj] / this->spacing_grid[1]  *  phi[2][ii] * gridpoint_potential[idx_z][idx_y][idx_x]; 
                                  force[2] +=  phi[0][kk] *  phi[1][jj] * dphi[2][ii] /  this->spacing_grid[2] * gridpoint_potential[idx_z][idx_y][idx_x];

                               }
                           }
                       }
                       FloatType tmp = config->sigma3[particle_type[particle_idx]] * config->sqrt_epsilon[particle_type[particle_idx]];

                       //update particle force
                       particle_force_grid[particle_idx][0] -= 4.0 * force[0] * tmp;
                       particle_force_grid[particle_idx][1] -= 4.0 * force[1] * tmp;
                       particle_force_grid[particle_idx][2] -= 4.0 * force[2] * tmp;

                       potential -= 2.0 * tmp * tmp_potential;
                    }// end particle_idx
                }
            }
        }
        for(int ii = 0; ii < 3; ++ii){
           delete[] phi[ii];
           delete[] dphi[ii];
        }
        return potential;
    }


    /*-----------------------------------------------------------------------------
     *  Class: CoarsestGrid
     *-----------------------------------------------------------------------------*/

    void msm::CoarsestGrid::init(msm::Config* config, int level)
    {
        this->spacing_grid[0] = config->grid_spacing[0] * pow(2.0,level);
        this->spacing_grid[1] = config->grid_spacing[1] * pow(2.0,level);
        this->spacing_grid[2] = config->grid_spacing[2] * pow(2.0,level);

        assert(config->num_gridpoints[0] % (1<<level) == 0);
        assert(config->num_gridpoints[1] % (1<<level) == 0);
        assert(config->num_gridpoints[2] % (1<<level) == 0);
        this->num_gridpoints[0] = config->num_gridpoints[0]/(1<<level);
        this->num_gridpoints[1] = config->num_gridpoints[1]/(1<<level);
        this->num_gridpoints[2] = config->num_gridpoints[2]/(1<<level);

        this->interpolation_order = config->interpolation_order;
        this->level               = level;

        allocate3DArray<FloatType> (this->gridpoint_potential, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        allocate3DArray<FloatType> (this->gridpoint_force_x, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        allocate3DArray<FloatType> (this->gridpoint_force_y, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        allocate3DArray<FloatType> (this->gridpoint_force_z, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        allocate3DArray<FloatType> (this->gridpoint_xi, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        allocate3DArray<FloatType> (this->direct_stencil, num_gridpoints[0], num_gridpoints[1], num_gridpoints[2]);
        
        // for coarsest grid not necessary; it would always be 1 * num_gridpoints[]
        this->direct_periodic_int_mult_of_num_gridpoints[0] = -1;
        this->direct_periodic_int_mult_of_num_gridpoints[1] = -1;
        this->direct_periodic_int_mult_of_num_gridpoints[2] = -1;

#ifdef PROFILE
        double start_time = omp_get_wtime();
#endif
        FloatType error_tol = 1.e-4;
        //parameters where iteration starts and how big the stepsize is
        //need to be clever choices!
        int maxiter = 8; //has to be larger or equal to two
        assert(maxiter >= 2);
        int iter_init = 50;
        int iterstepsize = 20;
        //Loop over all grid points
        printf("STARTING PRECALCULATION\n");
#pragma omp parallel for schedule(guided)
        for(int i=0; i < this->num_gridpoints[2]; ++i){
           FloatType dist_z = i * this->spacing_grid[2];
           for(int j=0;j < this->num_gridpoints[1]; ++j){
              FloatType dist_y = j * this->spacing_grid[1];
              for(int k=0;k < this->num_gridpoints[0]; ++k){
                 FloatType dist_x = k * this->spacing_grid[0];

                 FloatType old_value = 0.0;
                 this->direct_stencil[i][j][k] = this->seriescalc(dist_x,dist_y,dist_z, config->domain_extend, iter_init, config->cutoff, config->split_order);

                 for(int iteration = 1; iteration < maxiter; ++iteration){
                    int maxitseries = iter_init + iteration * iterstepsize;
                    old_value = this->direct_stencil[i][j][k];
                    this->direct_stencil[i][j][k] = this->seriescalc(dist_x, dist_y, dist_z, config->domain_extend, maxitseries, config->cutoff, config->split_order);
                    if(fabs((old_value - this->direct_stencil[i][j][k])/this->direct_stencil[i][j][k]) <= error_tol){
#ifdef SERIES_OUTPUT
                          printf("series converged for indices %d %d %d with %d images in each direction \n", i, j ,k , maxitseries );
#endif
                       break;
                    }
                 }
                 
#ifdef SERIES_OUTPUT
                   if(fabs((old_value - this->direct_stencil[i][j][k])/this->direct_stencil[i][j][k]) > error_tol)
                    printf("Not converged for indices %d %d %d (%e)\n", i, j ,k,
                          fabs((old_value - this->direct_stencil[i][j][k])/this->direct_stencil[i][j][k]));
#endif
              }
           }
        }
#ifdef PROFILE
        printf("Time spent in CoarsestGrid::init() %f sec\n\n",omp_get_wtime() - start_time);
#endif

    }

FloatType msm::CoarsestGrid::seriescalc(FloatType dist_x,FloatType dist_y,FloatType dist_z, FloatType* domain_extend, int max_iterations, FloatType cutoff, int split_order)
{
   FloatType seriesvalue = 0.;
   FloatType a2grid = pow(2.,2 * level) * sqr(cutoff);

   for(int i = (2 * max_iterations + 1); i > 0; --i){
      int sign_z = (i - (i/2) * 2) * 2 - 1;
      int image_ind_z = (i/2) * sign_z;
      FloatType image_dist_z = dist_z + image_ind_z * domain_extend[2]; 
      for(int j = (2 * max_iterations + 1); j > 0; --j){
         int sign_y = (j - (j/2) * 2) * 2 - 1;
         int image_ind_y = (j/2) * sign_y;
         FloatType image_dist_y = dist_y + image_ind_y * domain_extend[1];
         for(int k= (2 * max_iterations + 1); k > 0; --k){
            int sign_x = (k - (k/2) * 2) * 2 - 1;
            int image_ind_x = (k/2) * sign_x;
            FloatType image_dist_x = dist_x + image_ind_x * domain_extend[0];
            FloatType distance2 = image_dist_x * image_dist_x + image_dist_y * image_dist_y + image_dist_z * image_dist_z;
            FloatType distscaled2 = distance2 / a2grid; 
            if(distscaled2 <= 1.)
               seriesvalue += msm::functions::gamma(distscaled2, split_order);
            else
               seriesvalue += 1. / (distscaled2 * distscaled2 * distscaled2) ;
         }
      }
   }
   return seriesvalue / (a2grid*a2grid*a2grid);
}

    void msm::CoarsestGrid::direct(FloatType cutoff, FloatType *gamma_grid, int **gamma_grid_index, int gamma_grid_size){

        assert(this->gridpoint_potential[0][0][0] == 0);

        //iterate over all gridpoints
        for(int i=0;i < this->num_gridpoints[2] ;++i){
            for(int j=0;j < this->num_gridpoints[1] ;++j){
                for(int k=0;k < this->num_gridpoints[0] ;++k){

                   for(int ii = 0; ii < this->num_gridpoints[2]; ++ii){ 
                      int ii_tmp = (ii + i + this->num_gridpoints[2])%this->num_gridpoints[2]; 
                      for(int jj = 0; jj < this->num_gridpoints[1]; ++jj){ 
                         int jj_tmp = (jj + j + this->num_gridpoints[1])%this->num_gridpoints[1];
                         for(int kk = 0; kk < this->num_gridpoints[0]; ++kk){ 
                            //periodic boundaries
                            int kk_tmp = (kk + k + this->num_gridpoints[0])%this->num_gridpoints[0];

                            this->gridpoint_potential[i][j][k] += this->direct_stencil[ii][jj][kk] * this->gridpoint_xi[ii_tmp][jj_tmp][kk_tmp];
                         }
                      }
                   }
                }
            }
        }
    }
    
/*-----------------------------------------------------------------------------
 *  namespace: ShortRange
 *-----------------------------------------------------------------------------*/
    /**
     * \param[in] input_file_name String containing the file_name which contains
     *            the initial positions,forces, old forces, velocity and type of every
     *            particle 
     * \param[in] num_cells number of cells in each dimension. 0=x, 1=y, 2=z.
     * \param[in] domain_extend extend of the domain in each dimension.
     * \param[in] cutoff cutoff radius
     * \param[in] dt timestep for the time integration
     **/
    msm::ShortRange::ShortRange(
                 std::string   input_filename,
                 msm::Config* config)
    {
        this->config = config;
        this->potential_energy = 0.0;

        this->tmp_particle_id = new int[config->num_particles];
        this->tmp_particle_type = new int[config->num_particles];
        memory::create(this->tmp_particle_pos, config->num_particles, 3);
        memory::create(this->tmp_particle_force, config->num_particles, 3);
        memory::create(this->tmp_particle_force_old, config->num_particles, 3);
        memory::create(this->tmp_particle_vel, config->num_particles, 3);

        this->num_cells[0] = config->num_cells[0];
        this->num_cells[1] = config->num_cells[1];
        this->num_cells[2] = config->num_cells[2];

        this->read_dump(input_filename, config);
#ifdef DEBUG
        this->test();
#endif

        this->num_nbr_cells[0] = ceil(config->cutoff/this->cell_extend[0]);
        this->num_nbr_cells[1] = ceil(config->cutoff/this->cell_extend[1]);
        this->num_nbr_cells[2] = ceil(config->cutoff/this->cell_extend[2]);
        if(!(this->num_nbr_cells[0] == 1 && this->num_nbr_cells[1] == 1 && this->num_nbr_cells[2] == 1))
        {
            std::cout << "the cutoff became larger than the cell width in at least one dimension. Increase the cell width or decrease the cutoff." << std::endl;
            exit(-1);
        }
        this->computeSelfEnergy(config);
    }

    msm::ShortRange::~ShortRange()
    {
       delete[] tmp_particle_id;
       delete[] tmp_particle_type;
       delete[] &(tmp_particle_pos[0][0]);
       delete[] &(tmp_particle_force[0][0]);
       delete[] &(tmp_particle_force_old[0][0]);
       delete[] &(tmp_particle_vel[0][0]);
       free3DArray(this->cells, this->num_cells[1], this->num_cells[2]);
    }

void msm::ShortRange::computeSelfEnergy(msm::Config* config)
{
   this->self_energy = 0.0; //delete old value or initializing to zero
   //compute the self interaction energy of all particles without the constant term
   for(int i=0;i < num_cells[2]; ++i){
      for(int j=0;j < num_cells[1]; ++j){
         for(int k=0;k < num_cells[0]; ++k){
            this->self_energy += this->cells[i][j][k].computeSelfEnergy(config);
         }
      }
   }
   this->self_energy *= 2. * msm::functions::gamma(0., config->split_order) / ( sqr(sqr(config->cutoff)) * sqr(config->cutoff) ); //correction with the constant term
}
void msm::ShortRange::print()
{
   std::ofstream f;
   f.open ("output_initial_forces.txt");
   if( f == NULL ) { printf("WARNING: Can not open file!\n"); exit(-1); }
   f.precision(16);
   f << "id \t fx \t fy \t fz\n";
   //std::cout.precision(16);
   //std::cout << "id \t fx_total \t fy_total \t fz_total \t fx_cell \t fy_cell \t fz_cell \t fx_grid \t fy_grid \t fz_grid\n";

   for(int i=0;i < num_cells[2]; ++i)
      for(int j=0;j < num_cells[1]; ++j)
         for(int k=0;k < num_cells[0]; ++k)
            this->cells[i][j][k].print(&f);

   f.close();
}

/**
 * Reads data from a lammps generated dump file.
 * If particles are outside the domain, they will be remapped into the domain.
 * Moreover, the domain origin is shifted to (0,0,0).
 * The dumpfile is expected to have the columns in order: ID, particle type, particle mass, X Position, Y Position, Z Position. 
 */
void msm::ShortRange::read_dump(std::string input_filename, msm::Config* config)
{
   using namespace std;
   vector<string> tmp_strings;
   ifstream inputfile(input_filename.c_str());

   FloatType boxlo[3]; 
   FloatType boxhi[3];

   FloatType force[3]     = {0.0, 0.0, 0.0};
   FloatType force_old[3] = {0.0, 0.0, 0.0};
   FloatType vel[3]       = {0.0, 0.0, 0.0};
   FloatType pos[3]       = {0.0, 0.0, 0.0};

   int tmp_num_particles = 3;
   if(inputfile)
   {
      int i;
      string line;
      getline(inputfile, line);
      if(line.compare("ITEM: TIMESTEP") != 0 )
      {
         cerr << "incompatible input file!!! " << input_filename << " expected  ITEM: TIMESTEP in first line" << endl;
         exit(-1);
      }
      getline(inputfile, line);
      cout << "reading timestep # " << line << " (not further used for this simulation)"<< endl;
      getline(inputfile, line);
      if(line.compare("ITEM: NUMBER OF ATOMS") != 0 )
      {
         cerr << "incompatible input file!!! " << input_filename << " expected  ITEM: NUMBER OF ATOMS  in third line" << endl;
         exit(-1);
      }
      inputfile >> tmp_num_particles;
      assert(tmp_num_particles == this->config->num_particles);
      if(!inputfile)
      {
         cerr << "Could not read number of atoms from file!"<< endl;
         exit(-1);
      }
      // get rid of whitespace
      inputfile.ignore(1,' ');//std::numeric_limits<std::streamsize>::max(),' ');

      //allow each cell to store up to 10 times as mean particles as an average cell
      int max_particles_per_cell = this->config->num_particles / (this->num_cells[0] * this->num_cells[1] * this->num_cells[2]) * 10;
      max_particles_per_cell = max(max_particles_per_cell, 16); // each cell can store at least 16 particles

      allocate3DArray<Cell>(this->cells, this->num_cells[0], this->num_cells[1], this->num_cells[2]);
      for(int i=0;i < num_cells[2]; ++i){
         for(int j=0;j < num_cells[1]; ++j){
            for(int k=0;k < num_cells[0]; ++k){ //x-dim
               this->cells[i][j][k].init(max_particles_per_cell);
            }
         }
      }

      getline(inputfile, line);
      if(line.compare("ITEM: BOX BOUNDS pp pp pp") != 0 )
      {
         cerr << "incompatible input file!!! " << input_filename << " expected BOX BOUNDS pp pp pp in fifth line" << endl;
         exit(-1);
      }
      cout << "using periodic boundary conditions in each direction" << endl;
      for(i = 0; i < 3; i++ )
         inputfile >> boxlo[i] >> boxhi[i];
      if(!inputfile)
      {
         cerr << "Could not read domain extend from file!"<< endl;
         exit(-1);
      }
      cout << "box dimensions are:" << endl; 
      cout  << "domain: " << endl;
      FloatType tmp_domain_extend[3];
      for(i = 0; i < 3; i++ ){
         tmp_domain_extend[i] = boxhi[i] - boxlo[i];
         assert(config->domain_extend[i] == tmp_domain_extend[i]);
         cout  << config->domain_extend[i] <<" ";
      }
      cout  << "\ncell extend: " << endl;
      for(i = 0; i < 3; i++ ){
         this->cell_extend[i] = (this->config->domain_extend[i] / num_cells[i]);
         cout  <<this->cell_extend[i] <<" ";
      }
      cout  << endl;

      // get rid of whitespace
      inputfile.ignore(1,' ');


      getline(inputfile, line);
      //tmp_strings = split(line);
      //for(std::vector<string>::iterator iter = tmp_strings.begin(); iter != tmp_strings.end(); iter++)
      //  cout << *iter << endl;
      line.erase(find_if(line.rbegin(), line.rend(), not1(ptr_fun<int, int>(isspace))).base(), line.end()); // get rid of whitespace at the end of line

      for(i=0; i < this->config->num_particles; i++)
      {
         int id, type;
         if(line.compare("ITEM: ATOMS id type x y z vx vy vz fx fy fz") == 0 ){
             inputfile >> id >> type >> pos[0] >> pos[1] >> pos[2] >> vel[0] >> vel[1] >> vel[2] >> force_old[0] >> force_old[1] >> force_old[2];
         }else if(line.compare("ITEM: ATOMS id type x y z") == 0 ){
             inputfile >> id >> type >> pos[0] >> pos[1] >> pos[2];
         }else{
             cerr << "incompatible data format of input file " << input_filename << "!!!" << endl << "id type x y z vx vy vz fx fy fz" << endl << "OR id type x y z" << endl << "expected!" << endl;
             exit(-1);
         }
         //transform type from starting at 1 to start at 0
         --type;
         int j;
         // shift particle positions such that the domain origin is at (0,0,0)
         for(j=0;j < 3; ++j){
            pos[j] -= boxlo[j];
            //ensure that all particles are within the domain.
            while(pos[j] < 0.0)
               pos[j] += this->config->domain_extend[j];
            while(pos[j] > this->config->domain_extend[j])
               pos[j] -= this->config->domain_extend[j];
         }

         this->mapParticleToCell(id, type, force, force_old, vel, pos); //force_old = force
         if(!inputfile)
         {
            cerr << "Could not read the " << i << "-th atom of the inputfile!"<< endl;
            exit(-1);
         }
      }
   }
   else
   {
      cerr << "cannot open file " << input_filename << endl;
      exit(-1);
   }
}


    /**
     *  For geometricHashing in the course of the simulation please refer to updateGeometricHashing()
     **/
    void msm::ShortRange::mapParticleToCell( 
          int id,
          int type,
          FloatType force[3], 
          FloatType force_old[3],
          FloatType vel[3],
          FloatType pos[3])
    {
        // compute index depending on pos
        int i = static_cast<int>( floor(pos[2]/this->cell_extend[2]) ); //z
        int j = static_cast<int>( floor(pos[1]/this->cell_extend[1]) ); //y
        int k = static_cast<int>( floor(pos[0]/this->cell_extend[0]) ); //x

        assert(i>= 0 && i < this->num_cells[2] 
              && j>= 0 && j < this->num_cells[1] 
              && k>= 0 && k < this->num_cells[0]) ;

        // update respective cell...
        Cell& responsible_cell = this->cells[i][j][k];
//        printf("cell:%d %d %d : %d\n",i,j,k,cells[i][j][k].getNumParticles());
       
        responsible_cell.addParticle(id, type, force, force_old, vel, pos);
    }

FloatType msm::ShortRange::kineticEnergy()
{
   FloatType energy = 0.0;
   //loop over all cells
   for(int i=0;i < this->num_cells[2] ;++i){
      for(int j=0;j < this->num_cells[1] ;++j){
         for(int k=0;k < this->num_cells[0] ;++k){
            Cell& current_cell = this->cells[i][j][k];
            energy += current_cell.kineticEnergy(this->config->particle_mass);
         }
      }
   }
   return energy;
}

    void msm::ShortRange::computeForcesInteriorNonSym()
{
   this->potential_energy = 0.0;

   //loop over all cells
   for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // y-dim
         for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               assert(i + ii >= 0 && i + ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  assert(j + jj >= 0 && j + jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                        assert(k + kk >= 0 && k + kk < this->num_cells[0]);
                        Cell& nbr_cell = this->cells[i + ii][j + jj][k + kk];
                        this->potential_energy += current_cell.computeInteractionWithCellNonSym( 
                              nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma, this->config->split_order );
                  }
               }
            } // current_cell done
         }
      }
   }
}


void msm::ShortRange::computeForcesInteriorSym()
{
   this->potential_energy = 0.0;

   //loop over all cells
   for(int i=num_nbr_cells[2]; i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
      for(int j=num_nbr_cells[1]; j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // y-dim
         for(int k=num_nbr_cells[0]; k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim

            Cell& current_cell = this->cells[i][j][k];
            int ii,jj,kk;
            ii = 1; 
            for(jj = -num_nbr_cells[1] ; jj < num_nbr_cells[1]; ++jj){
               for(kk = -num_nbr_cells[0] ; kk < num_nbr_cells[0]; ++kk){
                  Cell& nbr_cell = this->cells[i + ii][j + jj][k + kk];
                  this->potential_energy += current_cell.computeInteractionWithCellSym( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma, this->config->split_order );
               }
            }
            ii = 0; 
            jj = 1; 
            for(kk = -num_nbr_cells[0] ; kk < num_nbr_cells[0]; ++kk){
               Cell& nbr_cell = this->cells[i + ii][j + jj][k + kk];
               this->potential_energy += current_cell.computeInteractionWithCellSym( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma, this->config->split_order );
            }
            ii = 0; 
            jj = 0; 
            kk = 1;
            Cell& nbr_cell = this->cells[i + ii][j + jj][k + kk];
            this->potential_energy += current_cell.computeInteractionWithCellSym( 
                  nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma, this->config->split_order );

            this->potential_energy += current_cell.computeInteractionWithItselfSym(
                  this->config->cutoff_sqr, this->config->epsilon, this->config->sigma, this->config->split_order );
         }
      }
   }
}

void msm::ShortRange::computeForcesBoundary()
{
   int shift_coords[] = {0, 0, 0};

   /************************** 
    * Faces 
    **************************/ 

   //front
   {
      int i = 0;
      int tmp_ii;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // y-dim
         for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               //periodic boundarys
               if(ii < 0 ){
                  shift_coords[2] = -1;
                  tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
               }else{
                  shift_coords[2] = 0;
                  tmp_ii = i + ii;
               }
               assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  assert(j + jj >= 0 && j + jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                     assert(k + kk >= 0 && k + kk < this->num_cells[0]);
                     Cell& nbr_cell = this->cells[tmp_ii ][j + jj][k + kk];
                     this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                           nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
                  }
               }
            } // current_cell done
         }
      }
   }

   // back
   {
      int i = this->num_cells[2] - 1;
      int tmp_ii;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // y-dim
         for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               //periodic boundarys
               if(ii > 0 ){
                  shift_coords[2] = 1;
                  tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
               }else{
                  shift_coords[2] = 0;
                  tmp_ii = i + ii;
               }
               assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  assert(j + jj >= 0 && j + jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                     assert(k + kk >= 0 && k + kk < this->num_cells[0]);
                     Cell& nbr_cell = this->cells[tmp_ii ][j + jj][k + kk];
                     this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                           nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
                  }
               }
            } // current_cell done
         }
      }
   }

   //left 
   {
      int j = 0;
      int tmp_jj;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               assert(i + ii>= 0 && i + ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  //periodic boundarys
                  if(jj < 0 ){
                     shift_coords[1] = -1;
                     tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
                  }else{
                     shift_coords[1] = 0;
                     tmp_jj = j + jj;
                  }
                  assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                     assert(k + kk>= 0 && k + kk < this->num_cells[0]);
                     Cell& nbr_cell = this->cells[i + ii ][tmp_jj][k + kk];
                     this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                           nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
                  }
               }
            } // current_cell done
         }
      }
   }

   //right
   {
      int j = this->num_cells[1] - 1;
      int tmp_jj;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               assert(i + ii>= 0 && i + ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  //periodic boundarys
                  if(jj > 0 ){
                     shift_coords[1] = 1;
                     tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
                  }else{
                     shift_coords[1] = 0;
                     tmp_jj = j + jj;
                  }
                  assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                     assert(k + kk>= 0 && k + kk < this->num_cells[0]);
                     Cell& nbr_cell = this->cells[i + ii ][tmp_jj][k + kk];
                     this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                           nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
                  }
               }
            } // current_cell done
         }
      }
   }

   //bottom
   {
      int k = 0;
      int tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // y-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               assert(i + ii>= 0 && i + ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  assert(j + jj>= 0 && j + jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                     //periodic boundarys
                     if(kk < 0 ){
                        shift_coords[0] = -1;
                        tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                     }else{
                        shift_coords[0] = 0;
                        tmp_kk = k + kk;
                     }
                     assert(tmp_kk >= 0 &&tmp_kk < this->num_cells[0]);

                     Cell& nbr_cell = this->cells[i + ii ][j + jj][tmp_kk];
                     this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                           nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
                  }
               }
            } // current_cell done
         }
      }
   }

   //top
   {
      int k = this->num_cells[0] - 1;
      int tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // y-dim

            Cell& current_cell = this->cells[i][j][k];

            //loop over all neighboring cells 
            for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
               assert(i + ii>= 0 && i + ii < this->num_cells[2]);
               for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
                  assert(j + jj>= 0 && j + jj < this->num_cells[1]);
                  for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                     //periodic boundarys
                     if(kk > 0 ){
                        shift_coords[0] = 1;
                        tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                     }else{
                        shift_coords[0] = 0;
                        tmp_kk = k + kk;
                     }
                     assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);

                     Cell& nbr_cell = this->cells[i + ii ][j + jj][tmp_kk];
                     this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                           nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
                  }
               }
            } // current_cell done
         }
      }
   }
   
   /************************** 
    * edges
    **************************/ 
   // front left
   {
      int i = 0;
      int j = 0;
      int tmp_ii, tmp_jj;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii < 0 ){
               shift_coords[2] = -1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj < 0 ){
                  shift_coords[1] = -1;
                  tmp_jj = ((j+ jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j+ jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  assert(k + kk >= 0 && k + kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][k + kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // front right 
   {
      int i = 0;
      int j = num_cells[1] - 1;
      int tmp_ii, tmp_jj;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii < 0 ){
               shift_coords[2] = -1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj > 0 ){
                  shift_coords[1] = 1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  assert(kk + k >= 0 && kk + k < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][k + kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // front bottom 
   {
      int i = 0;
      int k = 0;
      int tmp_ii, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii < 0 ){
               shift_coords[2] = -1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               assert(j + jj >= 0 && j + jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk < 0 ){
                     shift_coords[0] = -1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k + kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][j + jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // front top
   {
      int i = 0;
      int k = num_cells[0] - 1;
      int tmp_ii, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii < 0 ){
               shift_coords[2] = -1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               assert(j + jj >= 0 && j + jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk > 0 ){
                     shift_coords[0] = 1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k + kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][j + jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // back left 
   {
      int i = num_cells[2] - 1;
      int j = 0;
      int tmp_ii, tmp_jj;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii > 0 ){
               shift_coords[2] = 1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <=num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj < 0 ){
                  shift_coords[1] = -1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  assert(k + kk >= 0 && k + kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][k + kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // back right 
   {
      int i = num_cells[2] - 1;
      int j = num_cells[1] - 1;
      int tmp_ii, tmp_jj;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int k=num_nbr_cells[0];k < this->num_cells[0] - num_nbr_cells[0] ;++k){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii > 0 ){
               shift_coords[2] = 1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj > 0 ){
                  shift_coords[1] = 1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  assert(k + kk >= 0 && k + kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][k + kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // back bottom 
   {
      int i = num_cells[2] - 1;
      int k = 0;
      int tmp_ii, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii > 0 ){
               shift_coords[2] = 1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               assert(j + jj>= 0 && j + jj< this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk < 0 ){
                     shift_coords[0] = -1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k + kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][j + jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // back top
   {
      int i = num_cells[2] - 1;
      int k = num_cells[0] - 1;
      int tmp_ii, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int j=num_nbr_cells[1];j < this->num_cells[1] - num_nbr_cells[1] ;++j){ // x-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            //periodic boundarys
            if(ii > 0 ){
               shift_coords[2] = 1;
               tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
            }else{
               shift_coords[2] = 0;
               tmp_ii = i + ii;
            }
            assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               assert(j + jj >= 0 && j + jj < this->num_cells[1]);
               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk > 0 ){
                     shift_coords[0] = 1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k + kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
                  Cell& nbr_cell = this->cells[tmp_ii ][j + jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // bottom left 
   {
      int j = 0;
      int k = 0;
      int tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            assert(i + ii >= 0 && i + ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj < 0 ){
                  shift_coords[1] = -1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);

               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk < 0 ){
                     shift_coords[0] = -1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k +kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);

                  Cell& nbr_cell = this->cells[i + ii ][tmp_jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // bottom right  
   {
      int j = num_cells[1] - 1;
      int k = 0;
      int tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            assert(i + ii >= 0 && i + ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj > 0 ){
                  shift_coords[1] = 1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);

               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk < 0 ){
                     shift_coords[0] = -1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k + kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);

                  Cell& nbr_cell = this->cells[i + ii ][tmp_jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // top left 
   {
      int j = 0;
      int k = num_cells[0] - 1;
      int tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            assert(i + ii >= 0 && i + ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj < 0 ){
                  shift_coords[1] = -1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);

               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk > 0 ){
                     shift_coords[0] = 1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k +kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);

                  Cell& nbr_cell = this->cells[i + ii ][tmp_jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   // top right  
   {
      int j = num_cells[1] - 1;
      int k = num_cells[0] - 1;
      int tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      for(int i=num_nbr_cells[2];i < this->num_cells[2] - num_nbr_cells[2] ;++i){ // z-dim
         Cell& current_cell = this->cells[i][j][k];
         //loop over all neighboring cells 
         for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
            assert(i + ii >= 0 && i + ii < this->num_cells[2]);
            for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
               //periodic boundarys
               if(jj > 0 ){
                  shift_coords[1] = 1;
                  tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
               }else{
                  shift_coords[1] = 0;
                  tmp_jj = j + jj;
               }
               assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);

               for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
                  //periodic boundarys
                  if(kk > 0 ){
                     shift_coords[0] = 1;
                     tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
                  }else{
                     shift_coords[0] = 0;
                     tmp_kk = k + kk;
                  }
                  assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);

                  Cell& nbr_cell = this->cells[i + ii ][tmp_jj][tmp_kk];
                  this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                        nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
               }
            }
         } // current_cell done
      }
   }
   /************************** 
    * Corners 
    **************************/ 
   // front left bottom
   {
      int i = 0;
      int j = 0;
      int k = 0;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii < 0 ){
            shift_coords[2] = -1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj < 0 ){
               shift_coords[1] = -1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk < 0 ){
                  shift_coords[0] = -1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // front left top 
   {
      int i = 0;
      int j = 0;
      int k = num_cells[0] - 1;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii < 0 ){
            shift_coords[2] = -1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj < 0 ){
               shift_coords[1] = -1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk > 0 ){
                  shift_coords[0] = 1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // front right bottom
   {
      int i = 0;
      int j = num_cells[1] - 1;
      int k = 0;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii < 0 ){
            shift_coords[2] = -1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj > 0 ){
               shift_coords[1] = 1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk < 0 ){
                  shift_coords[0] = -1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // front right top 
   {
      int i = 0;
      int j = num_cells[1] - 1;
      int k = num_cells[0] - 1;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii < 0 ){
            shift_coords[2] = -1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj > 0 ){
               shift_coords[1] = 1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk > 0 ){
                  shift_coords[0] = 1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // back left bottom
   {
      int i = num_cells[2] - 1;
      int j = 0;
      int k = 0;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii > 0 ){
            shift_coords[2] = 1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj < 0 ){
               shift_coords[1] = -1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk < 0 ){
                  shift_coords[0] = -1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // back left top 
   {
      int i = num_cells[2] - 1;
      int j = 0;
      int k = num_cells[0] - 1;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii > 0 ){
            shift_coords[2] = 1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj < 0 ){
               shift_coords[1] = -1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk > 0 ){
                  shift_coords[0] = 1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // back right bottom
   {
      int i = num_cells[2] - 1;
      int j = num_cells[1] - 1;
      int k = 0;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii > 0 ){
            shift_coords[2] = 1;
            tmp_ii = ((i +ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj > 0 ){
               shift_coords[1] = 1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk < 0 ){
                  shift_coords[0] = -1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
   // back right top 
   {
      int i = num_cells[2] - 1;
      int j = num_cells[1] - 1;
      int k = num_cells[0] - 1;
      int tmp_ii, tmp_jj, tmp_kk;
      shift_coords[0] = 0;
      shift_coords[1] = 0;
      shift_coords[2] = 0;
      Cell& current_cell = this->cells[i][j][k];
      //loop over all neighboring cells 
      for(int ii = -num_nbr_cells[2] ; ii <= num_nbr_cells[2]; ++ii){
         //periodic boundarys
         if(ii > 0 ){
            shift_coords[2] = 1;
            tmp_ii = ((i + ii + num_cells[2]) % num_cells[2]);
         }else{
            shift_coords[2] = 0;
            tmp_ii = i + ii;
         }
         assert(tmp_ii >= 0 && tmp_ii < this->num_cells[2]);
         for(int jj = -num_nbr_cells[1] ; jj <= num_nbr_cells[1]; ++jj){
            //periodic boundarys
            if(jj > 0 ){
               shift_coords[1] = 1;
               tmp_jj = ((j + jj + num_cells[1]) % num_cells[1]);
            }else{
               shift_coords[1] = 0;
               tmp_jj = j + jj;
            }
            assert(tmp_jj >= 0 && tmp_jj < this->num_cells[1]);
            for(int kk = -num_nbr_cells[0] ; kk <= num_nbr_cells[0]; ++kk){
               //periodic boundarys
               if(kk > 0 ){
                  shift_coords[0] = 1;
                  tmp_kk = ((k + kk + num_cells[0]) % num_cells[0]);
               }else{
                  shift_coords[0] = 0;
                  tmp_kk = k + kk;
               }
               assert(tmp_kk >= 0 && tmp_kk < this->num_cells[0]);
               Cell& nbr_cell = this->cells[tmp_ii ][tmp_jj][tmp_kk];
               this->potential_energy += current_cell.computeInteractionWithCellNonSymPeriodic( 
                     nbr_cell, this->config->cutoff_sqr, this->config->epsilon, this->config->sigma , shift_coords, this->config->domain_extend, this->config->split_order);
            }
         }
      } // current_cell done
   }
}

/**
 * \param[in] symmetric iff symmetric == 1 then use a symmetric version. Otherwise use nonsymmetric (default).
 * \return The particle_force array of each cell will contain the force
 *         solely due to the short-range interaction of the particles
 **/
void msm::ShortRange::computeForces(int symmetric)
{
   if(symmetric == 1)
      this->computeForcesInteriorSym();
   else
      this->computeForcesInteriorNonSym();
   this->computeForcesBoundary();
}

    /**
     * \return The particle_pos array of each cell will be updated 
     **/
    void msm::ShortRange::integrate(FloatType dt)
    {
        for(int i=0;i < this->num_cells[2] ;++i){
            for(int j=0;j < this->num_cells[1] ;++j){
                for(int k=0;k < this->num_cells[0] ;++k){
                    Cell& current_cell = this->cells[i][j][k];
                    // update particle pos and velocity
                    solver::verlet(  dt,
                                     current_cell.getNumParticles(),
                                     current_cell.getParticleType(),
                                     this->config->particle_mass,
                                     current_cell.getParticleForceCell(),
                                     current_cell.getParticleForceOld(), 
                                     current_cell.getParticleVel(), 
                                     current_cell.getParticlePos());
                }
            }
        }
    }
            
    void msm::ShortRange::periodicBoundary()
{
       for(int i=0;i < this->num_cells[2] ;++i){
          for(int j=0;j < this->num_cells[1] ;++j){
             for(int k=0;k < this->num_cells[0] ;++k){

                Cell& current_cell = this->cells[i][j][k];
                FloatType** particle_pos = current_cell.getParticlePos();

                for( int l = 0; l < current_cell.getNumParticles(); ++l){
                   while(particle_pos[l][0] < 0.0)
                      particle_pos[l][0] += this->config->domain_extend[0];
                   while(particle_pos[l][0] > this->config->domain_extend[0])
                      particle_pos[l][0] -= this->config->domain_extend[0];
                   while(particle_pos[l][1] < 0.0)
                      particle_pos[l][1] += this->config->domain_extend[1];
                   while(particle_pos[l][1] > this->config->domain_extend[1])
                      particle_pos[l][1] -= this->config->domain_extend[1];
                   while(particle_pos[l][2] < 0.0)
                      particle_pos[l][2] += this->config->domain_extend[2];
                   while(particle_pos[l][2] > this->config->domain_extend[2])
                      particle_pos[l][2] -= this->config->domain_extend[2];
                }
             }
          }
       }
}

    void msm::ShortRange::updateGeometricHashing() 
    {
       // extract particle data from cells
       int counter = 0;
       for(int i=0;i < this->num_cells[2] ;++i){
          for(int j=0;j < this->num_cells[1] ;++j){
             for(int k=0;k < this->num_cells[0] ;++k){

                Cell& current_cell = this->cells[i][j][k];
                int* particle_id = current_cell.getParticleID();
                int* particle_type = current_cell.getParticleType();
                FloatType** particle_pos = current_cell.getParticlePos();
                FloatType** particle_vel = current_cell.getParticleVel();
                FloatType** particle_force = current_cell.getParticleForceCell();
                FloatType** particle_force_old = current_cell.getParticleForceOld();

                for( int l = 0; l < current_cell.getNumParticles(); ++l){
                   this->tmp_particle_id[counter] = particle_id[l];
                   this->tmp_particle_type[counter] = particle_type[l];
                   this->tmp_particle_pos[counter][0] = particle_pos[l][0];
                   this->tmp_particle_pos[counter][1] = particle_pos[l][1];
                   this->tmp_particle_pos[counter][2] = particle_pos[l][2];
                   this->tmp_particle_vel[counter][0] = particle_vel[l][0];
                   this->tmp_particle_vel[counter][1] = particle_vel[l][1];
                   this->tmp_particle_vel[counter][2] = particle_vel[l][2];
                   this->tmp_particle_force[counter][0] = particle_force[l][0];
                   this->tmp_particle_force[counter][1] = particle_force[l][1];
                   this->tmp_particle_force[counter][2] = particle_force[l][2];
                   this->tmp_particle_force_old[counter][0] = particle_force_old[l][0];
                   this->tmp_particle_force_old[counter][1] = particle_force_old[l][1];
                   this->tmp_particle_force_old[counter][2] = particle_force_old[l][2];
                   ++counter;
                }

                current_cell.resetParticleCount();
             }
          }
       }
       //map particles to cells
       for(int i=0;i < this->config->num_particles ; ++i)
          this->mapParticleToCell(
                this->tmp_particle_id[i], 
                this->tmp_particle_type[i], 
                this->tmp_particle_force[i], 
                this->tmp_particle_force_old[i], 
                this->tmp_particle_vel[i], 
                this->tmp_particle_pos[i]);
    }
           
    void msm::ShortRange::swapForcePtr()
    {
        for(int i=0;i < this->num_cells[2] ;++i){
            for(int j=0;j < this->num_cells[1] ;++j){
                for(int k=0;k < this->num_cells[0] ;++k){
                    this->cells[i][j][k].swapForcePtr();
                    this->cells[i][j][k].resetForces();
                }
            }
        }
    }
 
    void msm::ShortRange::combineForces()
    {
        for(int i=0;i < this->num_cells[2] ;++i){
            for(int j=0;j < this->num_cells[1] ;++j){
                for(int k=0;k < this->num_cells[0] ;++k){
                    this->cells[i][j][k].combineForces();
                }
            }
        }
    }


    /*-----------------------------------------------------------------------------
     *  Class: Cell
     *-----------------------------------------------------------------------------*/
    msm::Cell::Cell()
    {
        this->num_particles = 0;  
        this->max_particles = 0; 

        this->particle_vel = NULL;
        this->particle_pos = NULL;
        this->particle_force_from_cells = NULL;
        this->particle_force_from_grids = NULL;
        this->particle_force_old = NULL;
        this->particle_type       = NULL;    
        this->particle_id = NULL;    
    }

    msm::Cell::~Cell()
    {
        delete[] this->particle_id;
        delete[] this->particle_type;
        delete[] this->particle_pos[0]; 
        delete[] this->particle_vel[0];
        delete[] this->particle_force_from_cells[0]; 
        delete[] this->particle_force_from_grids[0]; 
        delete[] this->particle_force_old[0];
    }

FloatType msm::Cell::computeSelfEnergy(msm::Config* config)
{
   FloatType e_self = 0.; // self interaction energy of the cell without constants
   FloatType tmp_xi; // temporary variable for the xi value of particle particle_idx
   for(int particle_idx=0; particle_idx < num_particles; ++particle_idx){
       tmp_xi = config->sigma3[this->particle_type[particle_idx]] * config->sqrt_epsilon[this->particle_type[particle_idx]];
       e_self += sqr(tmp_xi);
   }
   return e_self;
}

void msm::Cell::print(std::ofstream *output)
{
   for(int i=0;i < this->num_particles ;++i){
      *output << this->particle_id[i] << "\t";
      *output << this->particle_force_from_cells[i][0] << "\t";
      *output << this->particle_force_from_cells[i][1] << "\t";
      *output << this->particle_force_from_cells[i][2] << std::endl;
      /*std::cout << this->particle_id[i] << "\t";
       std::cout << this->particle_force_from_cells[i][0]+this->particle_force_from_grids[i][0] << "\t";
       std::cout << this->particle_force_from_cells[i][1]+this->particle_force_from_grids[i][1] << "\t";
       std::cout << this->particle_force_from_cells[i][2]+this->particle_force_from_grids[i][2] << "\t";
       std::cout << this->particle_force_from_cells[i][0] << "\t";
       std::cout << this->particle_force_from_cells[i][1] << "\t";
       std::cout << this->particle_force_from_cells[i][2] << "\t";
       std::cout << this->particle_force_from_grids[i][0] << "\t";
       std::cout << this->particle_force_from_grids[i][1] << "\t";
       std::cout << this->particle_force_from_grids[i][2] << std::endl;*/
   }
}
    /**
     * \param domain_extend domain extend in x,y and z direction. Domain origin is always at (0,0,0)
     */
void msm::Cell::test(FloatType* particle_mass, FloatType domain_extend[3])
{
   for(int i=0;i < this->num_particles ;++i){
      if(particle_mass[this->particle_type[i]] <= 0.0){
         std::cerr<< "ERROR: Particle mass <= 0.0\n";
         exit(-1);
      }
      if(this->particle_type[i] < 0){
         std::cerr<< "ERROR: Particle has invalid type \n";
         exit(-1);
      }
      if(this->particle_pos[i][0] < 0.0 || this->particle_pos[i][1] < 0.0 || this->particle_pos[i][2] < 0.0
            ||this->particle_pos[i][0] > domain_extend[0] 
            || this->particle_pos[i][1] > domain_extend[1] 
            || this->particle_pos[i][2] > domain_extend[2]){
         std::cerr<< "ERROR: Particle has left the domain!\n";
         exit(-1);
      }
   }
#ifdef DEBUG_FULL
   this->test_kineticEnergy();
   exit(-1);
#endif
}

/**
 * WARNING: This function changes the values for the mass and velocity of each particle.
 * After calling this function make sure to exit the simulation
 */
void msm::Cell::test_kineticEnergy(FloatType* particle_mass, int num_particle_types )
{
   //initialize with random values
   for(int i=0;i < this->num_particles ;++i){
      this->particle_type[i] = rand()%this->num_particle_types;
      this->particle_vel[i][0] = rand()/((FloatType)RAND_MAX);
      this->particle_vel[i][1] = rand()/((FloatType)RAND_MAX);
      this->particle_vel[i][2] = rand()/((FloatType)RAND_MAX);
   }
   for(int i=0;i < num_particle_types ;++i){
      particle_mass[i] = rand()/((FloatType)RAND_MAX);
   }

   FloatType e1 = this->kineticEnergy(particle_mass);
   //particles that are twice as heavy have to result in twice the kinetic energy
   for(int i=0;i < this->num_particle_types ;++i)
      particle_mass[i] *= 2;
   FloatType e2 = this->kineticEnergy(particle_mass);
   if(e1*2.0 - e2 > 1e-15){
      std::cerr<< "ERROR: kineticEngery function seems to be buggy!\n";
      exit(-1);
   }
   //chaning velocity components may not change the result
   for(int i=0;i < this->num_particles ;++i){
      FloatType tmp;
      tmp = this->particle_vel[i][0];
      this->particle_vel[i][0] = this->particle_vel[i][1];
      this->particle_vel[i][1] = tmp;
   }
   FloatType e3 = this->kineticEnergy(particle_mass);
   if(e3 - e2 > 1e-15){
      std::cerr<< "ERROR: kineticEngery function seems to be buggy!\n";
      exit(-1);
   }
   //particles that are sqrt(8.0) as fast have to result in 8 times the kinetic energy
   for(int i=0;i < this->num_particles ;++i){
      this->particle_vel[i][0] *= 2.0;
      this->particle_vel[i][1] *= 2.0;
      this->particle_vel[i][2] *= 2.0;
   }
   FloatType e4 = this->kineticEnergy(particle_mass);
   if(e3 * 4.0 - e4 > 1e-15){
      std::cerr<< "ERROR: kineticEngery function seems to be buggy!\n";
      exit(-1);
   }
   delete[] particle_mass;
}
    void msm::Cell::init(int max_particles)
    {  
        this->max_particles = max_particles;

        //allocate data contiguous in memory
        memory::create(this->particle_pos, this->max_particles, 3);
        memory::create(this->particle_vel, this->max_particles, 3);
        memory::create(this->particle_force_from_cells, this->max_particles, 3);
        memory::create(this->particle_force_from_grids, this->max_particles, 3);
        memory::create(this->particle_force_old, this->max_particles, 3);
        this->particle_type             = new int[this->max_particles];
        this->particle_id = new int[this->max_particles];
    }

    /**
     * \param[in] force     force of the current time step in x,y and z-direction
     * \param[in] force_old force of the previous time stepin x,y and z-direction
     * \param[in] vel       velocity in x,y and z-direction
     * \param[in] pos       position in x,y and z-direction
     * */
    void msm::Cell::addParticle(  
          int id,
          int type,
          FloatType force[3], 
          FloatType force_old[3],
          FloatType vel[3],
          FloatType pos[3])
    {
        if(this->num_particles < this->max_particles){
            
            this->particle_force_from_cells[this->num_particles][0]  = force[0];
            this->particle_force_from_cells[this->num_particles][1]  = force[1];
            this->particle_force_from_cells[this->num_particles][2] = force[2];

            this->particle_force_from_grids[this->num_particles][0] = 0.0;
            this->particle_force_from_grids[this->num_particles][1] = 0.0;
            this->particle_force_from_grids[this->num_particles][2] = 0.0;

            this->particle_force_old[this->num_particles][0] = force_old[0];
            this->particle_force_old[this->num_particles][1] = force_old[1];
            this->particle_force_old[this->num_particles][2] = force_old[2];

            this->particle_vel[this->num_particles][0] = vel[0];
            this->particle_vel[this->num_particles][1] = vel[1];
            this->particle_vel[this->num_particles][2] = vel[2];

            this->particle_pos[this->num_particles][0] = pos[0];
            this->particle_pos[this->num_particles][1] = pos[1];
            this->particle_pos[this->num_particles][2] = pos[2];

            this->particle_type[this->num_particles] = type;
            this->particle_id[this->num_particles] = id;

            ++(this->num_particles);

        }else{
            std::cout << "ERROR: not enought space to fit all particles within this cell.\n";
            exit(-1); 
        }
    }

    /**
     * Updates the forces for each particle in this cell by accounting for all
     * the particle particle interaction with the provided neighboring cell. 
     *
     * \param[in] nbr_cell Neighboring cell
     * \param[in] cutoff cutoff radius after which the particle-particle
     *            interaction is zero
     * \return Potential Energy
     **/
    FloatType msm::Cell::computeInteractionWithCellNonSym(Cell& nbr_cell, const FloatType cutoff_sqr, FloatType** epsilon, FloatType** sigma, int split_order)
    {
        int* nbr_particle_type = nbr_cell.getParticleType();
        FloatType** nbr_particle_pos = nbr_cell.getParticlePos();
        FloatType cutoff_6 = sqr(cutoff_sqr) * cutoff_sqr;
        FloatType cutoff_8 = sqr(sqr(cutoff_sqr));
        FloatType epot = 0.0;

        for(int i=0;i < this->num_particles ;++i){
            for(int j=0;j < nbr_cell.getNumParticles() ;++j){

                FloatType dx = (nbr_particle_pos[j][0] - this->particle_pos[i][0]);
                FloatType dy = (nbr_particle_pos[j][1] - this->particle_pos[i][1]);
                FloatType dz = (nbr_particle_pos[j][2] - this->particle_pos[i][2]);
                FloatType distance_sqr = dx*dx + dy*dy + dz*dz;

                if (distance_sqr < cutoff_sqr && distance_sqr > 0.0) { //only compute interaction if particles are within cutoff 

                   FloatType epsilon_tmp = epsilon[this->particle_type[i]][nbr_particle_type[j]];
                   FloatType sigma_sqr = sigma[this->particle_type[i]][nbr_particle_type[j]];
                   sigma_sqr = sqr(sigma_sqr);
                   FloatType sigma_distance_sqr = sigma_sqr / distance_sqr;
                   FloatType sigma_distance_sqr_6 = sigma_distance_sqr * sigma_distance_sqr * sigma_distance_sqr ; // (sigma/distance)^6
                   FloatType sigma_distance_sqr_8 = sigma_distance_sqr * sigma_distance_sqr_6; // (sigma/distance)^8
                   FloatType sigma_distance_sqr_14 = sigma_distance_sqr_8 * sigma_distance_sqr_6; // (sigma/distance)^8
                   
                   epot += 2.0 * epsilon_tmp * (sqr(sigma_distance_sqr_6) - sigma_distance_sqr_6); // Lennard-Jones update
                   epot += 2.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_6 * msm::functions::gamma(distance_sqr/cutoff_sqr, split_order); // Multilevel correction

                   FloatType force = 24.0 * epsilon_tmp / sigma_sqr * (sigma_distance_sqr_8 - 2.0 * sigma_distance_sqr_14); 
                   force +=  4.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_8 * msm::functions::dgamma(distance_sqr/cutoff_sqr, split_order);
                   this->particle_force_from_cells[i][0] += force * dx;
                   this->particle_force_from_cells[i][1] += force * dy;
                   this->particle_force_from_cells[i][2] += force * dz;
                }
            }
        }
        return epot;
    }

    /**
     * Updates the forces for each particle in this cell by accounting for all
     * the particle particle interaction with the provided neighboring cell. 
     *
     * \param[in] nbr_cell Neighboring cell
     * \param[in] cutoff cutoff radius after which the particle-particle
     *            interaction is zero
     * \param[in] shift_coords shift_coords[0] == -1 iff x-coordinates need shifted by domain_extend[0] in negative x-dimention
     *                         shift_coords[0] == 1 iff x-coordinates need shifted by domain_extend[0] in positive x-dimention
     *                         shift_coords[1] == -1 iff y-coordinates need shifted by domain_extend[1] in negative y-dimention
     *                         shift_coords[1] == 1 iff y-coordinates need shifted by domain_extend[1] in positive y-dimention
     *                         shift_coords[2] == -1 iff z-coordinates need shifted by domain_extend[2] in negative z-dimention
     *                         shift_coords[2] == 1 iff z-coordinates need shifted by domain_extend[2] in positive z-dimention
     * \param[in] domain_extend Domain extend in x(idx:0),y(idx:1) and z(idx:2) direction 
     * \return Potential Energy
     **/
    FloatType msm::Cell::computeInteractionWithCellNonSymPeriodic(Cell& nbr_cell, const FloatType cutoff_sqr, FloatType** epsilon, FloatType** sigma, 
          int shift_coords[3], 
          FloatType domain_extend[3], int split_order)
    {
        int* nbr_particle_type = nbr_cell.getParticleType();
        FloatType** nbr_particle_pos = nbr_cell.getParticlePos();
        FloatType cutoff_6 = sqr(cutoff_sqr) * cutoff_sqr;
        FloatType cutoff_8 = sqr(sqr(cutoff_sqr));
        FloatType epot = 0.0;

        for(int i=0;i < this->num_particles ;++i){
            for(int j=0;j < nbr_cell.getNumParticles() ;++j){

                FloatType dx = (nbr_particle_pos[j][0] + shift_coords[0] * domain_extend[0] - this->particle_pos[i][0]);
                FloatType dy = (nbr_particle_pos[j][1] + shift_coords[1] * domain_extend[1] - this->particle_pos[i][1]);
                FloatType dz = (nbr_particle_pos[j][2] + shift_coords[2] * domain_extend[2] - this->particle_pos[i][2]);
                FloatType distance_sqr = dx*dx + dy*dy + dz*dz;

                if (distance_sqr < cutoff_sqr && distance_sqr > 0.0) { //only compute interaction if particles are within cutoff 
                // if (distance_sqr > cutoff_sqr) distance_sqr = cutoff_sqr; // alternative if
                   FloatType epsilon_tmp = epsilon[this->particle_type[i]][nbr_particle_type[j]];
                   FloatType sigma_sqr = sigma[this->particle_type[i]][nbr_particle_type[j]];
                   sigma_sqr = sqr(sigma_sqr);
                   FloatType sigma_distance_sqr = sigma_sqr / distance_sqr;
                   FloatType sigma_distance_sqr_6 = sigma_distance_sqr * sigma_distance_sqr * sigma_distance_sqr ; // (sigma/distance)^6
                   FloatType sigma_distance_sqr_8 = sigma_distance_sqr * sigma_distance_sqr_6; // (sigma/distance)^8
                   FloatType sigma_distance_sqr_14 = sigma_distance_sqr_8 * sigma_distance_sqr_6; // (sigma/distance)^8
                   
                   epot += 2.0 * epsilon_tmp * (sqr(sigma_distance_sqr_6) - sigma_distance_sqr_6); // Lennard-Jones update
                   epot += 2.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_6 * msm::functions::gamma(distance_sqr/cutoff_sqr, split_order); // Multilevel correction

                   FloatType force = 24.0 * epsilon_tmp / sigma_sqr * (sigma_distance_sqr_8 - 2.0 * sigma_distance_sqr_14); 
                   force +=  4.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_8 * msm::functions::dgamma(distance_sqr/cutoff_sqr, split_order);
                   this->particle_force_from_cells[i][0] += force * dx;
                   this->particle_force_from_cells[i][1] += force * dy;
                   this->particle_force_from_cells[i][2] += force * dz;
                }
            }
        }
        return epot;
    }


    FloatType msm::Cell::computeInteractionWithCellSym(Cell& nbr_cell, const FloatType cutoff_sqr, FloatType** epsilon, FloatType** sigma, int split_order)
    {
        assert(&(nbr_cell) != this);

        int* nbr_particle_type = nbr_cell.getParticleType();
        FloatType** nbr_particle_pos = nbr_cell.getParticlePos();
        FloatType cutoff_6 = sqr(cutoff_sqr) * cutoff_sqr;
        FloatType cutoff_8 = sqr(sqr(cutoff_sqr));
        FloatType epot = 0.0;

        for(int i=0;i < this->num_particles ;++i){
            for(int j=0;j < nbr_cell.getNumParticles() ;++j){

                FloatType dx = (nbr_particle_pos[j][0] - this->particle_pos[i][0]);
                FloatType dy = (nbr_particle_pos[j][1] - this->particle_pos[i][1]);
                FloatType dz = (nbr_particle_pos[j][2] - this->particle_pos[i][2]);
                FloatType distance_sqr = dx*dx + dy*dy + dz*dz;

                if (distance_sqr < cutoff_sqr) { //only compute interaction if particles are within cutoff 
                // if (distance_sqr > cutoff_sqr) distance_sqr = cutoff_sqr; // alternative if
                   FloatType epsilon_tmp = epsilon[this->particle_type[i]][nbr_particle_type[j]];
                   FloatType sigma_sqr = sigma[this->particle_type[i]][nbr_particle_type[j]];
                   sigma_sqr = sqr(sigma_sqr);
                   FloatType sigma_distance_sqr = sigma_sqr / distance_sqr;
                   FloatType sigma_distance_sqr_6 = sigma_distance_sqr * sigma_distance_sqr * sigma_distance_sqr ; // (sigma/distance)^6
                   FloatType sigma_distance_sqr_8 = sigma_distance_sqr * sigma_distance_sqr_6; // (sigma/distance)^8
                   FloatType sigma_distance_sqr_14 = sigma_distance_sqr_8 * sigma_distance_sqr_6; // (sigma/distance)^8
                   
                   epot += 4.0 * epsilon_tmp * (sqr(sigma_distance_sqr_6) - sigma_distance_sqr_6); // Lennard-Jones update
                   epot += 4.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_6 * msm::functions::gamma(distance_sqr/cutoff_sqr, split_order); // Multilevel correction

                   FloatType force = 24.0 * epsilon_tmp / sigma_sqr * (sigma_distance_sqr_8 - 2.0 * sigma_distance_sqr_14); 
                   force +=  4.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_8 * msm::functions::dgamma(distance_sqr/cutoff_sqr, split_order);
                   this->particle_force_from_cells[i][0] += force * dx;
                   this->particle_force_from_cells[i][1] += force * dy;
                   this->particle_force_from_cells[i][2] += force * dz;
                   this->particle_force_from_cells[j][0] -= force * dx;
                   this->particle_force_from_cells[j][1] -= force * dy;
                   this->particle_force_from_cells[j][2] -= force * dz;
                }
            }
        }
        return epot;
    }

    FloatType msm::Cell::computeInteractionWithItselfSym(const FloatType cutoff_sqr, FloatType** epsilon, FloatType** sigma, int split_order)
    {
        FloatType cutoff_6 = sqr(cutoff_sqr) * cutoff_sqr;
        FloatType cutoff_8 = sqr(sqr(cutoff_sqr));
        FloatType epot = 0.0;

        for(int i=0;i < this->num_particles ;++i){
            for(int j=i+1;j < this->num_particles ;++j){

                FloatType dx = (this->particle_pos[j][0] - this->particle_pos[i][0]);
                FloatType dy = (this->particle_pos[j][1] - this->particle_pos[i][1]);
                FloatType dz = (this->particle_pos[j][2] - this->particle_pos[i][2]);
                FloatType distance_sqr = dx*dx + dy*dy + dz*dz;

                if (distance_sqr > cutoff_sqr) { //only compute interaction if particles are within cutoff

                   FloatType epsilon_tmp = epsilon[this->particle_type[i]][this->particle_type[j]];
                   FloatType sigma_sqr = sigma[this->particle_type[i]][this->particle_type[j]];
                   sigma_sqr = sqr(sigma_sqr);
                   FloatType sigma_distance_sqr = sigma_sqr / distance_sqr;
                   FloatType sigma_distance_sqr_6 = sigma_distance_sqr * sigma_distance_sqr * sigma_distance_sqr ; // (sigma/distance)^6
                   FloatType sigma_distance_sqr_8 = sigma_distance_sqr * sigma_distance_sqr_6; // (sigma/distance)^8
                   FloatType sigma_distance_sqr_14 = sigma_distance_sqr_8 * sigma_distance_sqr_6; // (sigma/distance)^8
                   
                   epot += 4.0 * epsilon_tmp * (sqr(sigma_distance_sqr_6) - sigma_distance_sqr_6); // Lennard-Jones update
                   epot += 4.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_6 * msm::functions::gamma(distance_sqr/cutoff_sqr, split_order); // Multilevel correction

                   FloatType force = 24.0 * epsilon_tmp / sigma_sqr * (sigma_distance_sqr_8 - 2.0 * sigma_distance_sqr_14); 
                   force +=  4.0 * epsilon_tmp * sqr(sigma_sqr)*sigma_sqr / cutoff_8 * msm::functions::dgamma(distance_sqr/cutoff_sqr, split_order);
                   this->particle_force_from_cells[i][0] += force * dx;
                   this->particle_force_from_cells[i][1] += force * dy;
                   this->particle_force_from_cells[i][2] += force * dz;
                   this->particle_force_from_cells[j][0] -= force * dx;
                   this->particle_force_from_cells[j][1] -= force * dy;
                   this->particle_force_from_cells[j][2] -= force * dz;
                }
            }
        }
        return epot;
    }
    void msm::Cell::resetForces()
    {
        for(int i=0;i < this->max_particles ;++i){
            this->particle_force_from_cells[i][0] = 0.0;
            this->particle_force_from_cells[i][1] = 0.0;
            this->particle_force_from_cells[i][2] = 0.0;

            this->particle_force_from_grids[i][0] = 0.0;
            this->particle_force_from_grids[i][1] = 0.0;
            this->particle_force_from_grids[i][2] = 0.0;
        }
    }

    void msm::Cell::swapForcePtr()
    {
        FloatType **tmp;
        tmp = this->particle_force_from_cells;
        this->particle_force_from_cells = this->particle_force_old;
        this->particle_force_old = tmp;
    }

    void msm::Cell::combineForces()
    {
        for(int i=0;i < this->num_particles ;++i){
            this->particle_force_from_cells[i][0] += this->particle_force_from_grids[i][0];
            this->particle_force_from_cells[i][1] += this->particle_force_from_grids[i][1];
            this->particle_force_from_cells[i][2] += this->particle_force_from_grids[i][2];

        }
    }
    FloatType msm::Cell::kineticEnergy(FloatType* particle_mass)
    {
       FloatType energy = 0.0;
       int i,d; 
       for(i=0;i < this->num_particles ;++i){
          FloatType vel = 0.0;
          for(d=0;d < 3; ++d)
             vel+= sqr(this->particle_vel[i][d]);

          energy += .5 * particle_mass[this->particle_type[i]] * vel;
       }
       return energy;
    }

    /*-----------------------------------------------------------------------------
     *  namespace: MSM
     *-----------------------------------------------------------------------------*/

    /// check whether all forces are zero
    void check_forces_are_zero(int num_cells[3], msm::Cell*** cells)
    {
        for(int i=0;i < num_cells[2]; ++i){
            for(int j=0;j < num_cells[1]; ++j){
                for(int k=0;k < num_cells[0]; ++k){
                    msm::Cell&  current_cell  = cells[i][j][k];
                    FloatType** particle_force_cell = current_cell.getParticleForceCell(); 
                    FloatType** particle_force_grid = current_cell.getParticleForceGrid();
                    for(int ii=0; ii < current_cell.getNumParticles(); ++ii ){
                        for(int jj=0; jj < 3; ++jj) {
                            assert(particle_force_cell[ii][jj] == 0.);
                            assert(particle_force_grid[ii][jj] == 0.);
                        }
                    }
                }
            }
        }

    }

    /// propagates particle in time
    void update(FloatType dt, msm::ShortRange* short_range, msm::LongRange* long_range)
    {
        short_range->combineForces();

#ifndef NOFORCE
        short_range->print();
#endif

        short_range->integrate(dt);        //integrate particle in time

        short_range->swapForcePtr();     //swaps the force and force_old ptr of each cell. force_old can be overwritten in the next time step.

        short_range->periodicBoundary();

        short_range->updateGeometricHashing();
    }

    /// computes one timestep of the simulation. This function is only used by msm::run
    void step(FloatType dt, msm::ShortRange* short_range, msm::LongRange* long_range)
    {
        // check whether all forces are set to zero before starting their new computation
        check_forces_are_zero(short_range->getNumCells(), short_range->getCells());

#ifdef PROFILE
        double ticks0,ticks1;
        ticks0 = omp_get_wtime();
#endif

        short_range->computeForces(0);

#ifdef PROFILE
        ticks1 = omp_get_wtime();
        ticks_short_range += ticks1 - ticks0;
        //printf("%f, ",ticks1-ticks0);
#endif

#ifdef PROFILE
        ticks0 = omp_get_wtime();
#endif
        long_range->computeForces(short_range->getNumCells(), short_range->getCells()); //sR.computeForces() needs to be completed befor this can start (this dependecy can be removed by adding an extra for long-range forces within the cell class
#ifdef PROFILE
        ticks1 = omp_get_wtime();
        ticks_long_range += ticks1 - ticks0;
        //printf("%f, ",ticks1-ticks0);
        ticks0 = omp_get_wtime();
#endif
        update(dt, short_range, long_range);
#ifdef PROFILE
        ticks1 = omp_get_wtime();
        ticks_update += ticks1 - ticks0;
        //printf("%f \n ",ticks1-ticks0);
#endif

    }

    /**
     * \param config Encapsulates all parameters for the simulation
     * \param short_range This parameter needs to be initialized before calling this function.
     * \param long_range This parameter needs to be initialized before calling this function.
     */
    void msm::run(msm::Config* config, msm::ShortRange* short_range, msm::LongRange* long_range)
    {
          printf("STARTING SIMULATION\n");
          printf("Timestep, kin. Energy, pot. Energy, total Energy\n");
          FloatType kinEnergy = 0.0, potEnergy = 0.0;
          for(int i=0;i < config->num_iterations ;++i){
              if(i % config->output_frequency == 0){
                kinEnergy = short_range->kineticEnergy();
              }
              step(config->dt, short_range, long_range);

              if(i % config->output_frequency  == 0){
                potEnergy = short_range->getPotentialEnergy();
                potEnergy += short_range->getSelfEnergy();
                potEnergy += long_range->getPotentialEnergy();
                printf("%7d %.16e %.16e %.16e \n", i, kinEnergy/config->num_particles, potEnergy/config->num_particles, (kinEnergy + potEnergy)/config->num_particles);
              }
          }
#ifdef PROFILE
        for(int i=0;i < config->num_grids; ++i)
          printf("Time in grid %d: %f\n",i,direct_grid_time[i]);

        printf("Total time in local part: %f sec (avg: %f ms)\n",ticks_short_range,1000. * ticks_short_range/config->num_iterations);
        printf("Total time in grid-based part: %f sec (avg: %f ms)\n",ticks_long_range,1000. * ticks_long_range/config->num_iterations);
        printf("Total time in update: %f sec (avg: %f ms)\n",ticks_update,1000. * ticks_update/config->num_iterations);
        printf("%1.2f %.16e %.16e %.16e %.16e\n", config->cutoff, potEnergy/config->num_particles, short_range->getPotentialEnergy()/config->num_particles, long_range->getPotentialEnergy()/config->num_particles, short_range->getSelfEnergy()/config->num_particles);//for output
#endif
    }

/*-----------------------------------------------------------------------------
 *  msm::functions
 *-----------------------------------------------------------------------------*/
//    FloatType msm::functions::gamma(FloatType distance_sqr)
//    {
//       FloatType distance_12 = sqr(sqr(distance_sqr)) * sqr(distance_sqr);
//       return 15.0/8.0 - 5.0/4.0 * distance_12 + 3.0/8.0 * sqr(distance_12);
//    }

FloatType msm::functions::gamma(const FloatType distance_sqr, int split_order)
{
   assert(distance_sqr <= 1.);
   FloatType g = msm::functions::gcons[split_order][0];
   //g is the polynom
   FloatType distance12 = sqr(sqr(distance_sqr)) * sqr(distance_sqr);
   FloatType tmp = distance12;
   for (int n = 1; n <= split_order; n++) {
      g += msm::functions::gcons[split_order][n] * tmp;
      tmp *= distance12;
   }
   return g;
}

FloatType msm::functions::dgamma(const FloatType distance_sqr, int split_order)
{
   assert(distance_sqr <= 1.);
   FloatType distance10 = sqr(sqr(distance_sqr)) * distance_sqr;
   FloatType distance12 = distance10 * distance_sqr;

   FloatType dg = 0.;
   FloatType tmp = distance10;
   for (int n=0; n<split_order; n++) {
      dg += msm::functions::dgcons[split_order][n] * tmp;
      tmp *= distance12;
   }
   return dg;
}


    FloatType msm::functions::phi(FloatType distance, int order) //WARNING: distance should not be a reference since it will be overwritten
{
   double phi = 0.0;
   distance = (distance < 0.0) ? -distance : distance; 
   double distance2 = sqr(distance);

   if (order == 4) {
      if (distance <= 1) {
         phi = (1.0 - distance) * (1.0 + distance - 1.5*distance2);
      } else if (distance <= 2) {
         phi = -0.5*(distance - 1.0)*(2.0 - distance)*(2.0 - distance);
      } else {
         phi = 0.0;
      }

   } else if (order == 6) {
      if (distance <= 1) {
         phi = (1.0 - distance2)*(2.0 - distance)*(6.0 + 3.0*distance -
               5.0*distance2)/12.0;
      } else if (distance <= 2) {
         phi = -(distance - 1.0)*(2.0 - distance)*(3.0 - distance)*
            (4.0 + 9.0*distance - 5.0*distance2)/24.0;
      } else if (distance <= 3) {
         phi = (distance - 1.0)*(distance - 2.0)*(3.0 - distance)*
            (3.0 - distance)*(4.0 - distance)/24.0;
      } else {
         phi = 0.0;
      }

   } else if (order == 8) {
      if (distance <= 1) {
         phi = (1.0 - distance2)*(4.0 - distance2)*(3.0 - distance)*
            (12.0 + 4.0*distance - 7.0*distance2)/144.0;
      } else if (distance <= 2) {
         phi = -(distance2 - 1.0)*(2.0 - distance)*(3.0 - distance)*
            (4.0 - distance)*(10.0 + 12.0*distance - 7.0*distance2)/240.0;
      } else if (distance <= 3) {
         phi = (distance - 1.0)*(distance - 2.0)*(3.0 - distance)*(4.0 - distance)*
            (5.0 - distance)*(6.0 + 20.0*distance - 7.0*distance2)/720.0;
      } else if (distance <= 4) {
         phi = -(distance - 1.0)*(distance - 2.0)*(distance - 3.0)*(4.0 - distance)*
            (4.0 - distance)*(5.0 - distance)*(6.0 - distance)/720.0;
      } else {
         phi = 0.0;
      }

   } else if (order == 10) {
      if (distance <= 1) {
         phi = (1.0 - distance2)*(4.0 - distance2)*(9.0 - distance2)*
            (4.0 - distance)*(20.0 + 5.0*distance - 9.0*distance2)/2880.0;
      } else if (distance <= 2) {
         phi = -(distance2 - 1.0)*(4.0 - distance2)*(3.0 - distance)*(4.0 - distance)*
            (5.0 - distance)*(6.0 + 5.0*distance - 3.0*distance2)/1440.0;
      } else if (distance <= 3) {
         phi = (distance2 - 1.0)*(distance - 2.0)*(3.0 - distance)*(4.0 - distance)*
            (5.0 - distance)*(6.0 - distance)*(14.0 + 25.0*distance - 9.0*distance2)/10080.0;
      } else if (distance <= 4) {
         phi = -(distance - 1.0)*(distance - 2.0)*(distance - 3.0)*(4.0 - distance)*
            (5.0 - distance)*(6.0 - distance)*(7.0 - distance)*
            (8.0 + 35.0*distance - 9.0*distance2)/40320.0;
      } else if (distance <= 5) {
         phi = (distance - 1.0)*(distance - 2.0)*(distance - 3.0)*
            (distance - 4.0)*(5.0 - distance)*(5.0 - distance)*(6.0 - distance)*
            (7.0 - distance)*(8.0 - distance)/40320.0;
      } else {
         phi = 0.0;
      }
   }

   return phi;
}

    FloatType msm::functions::dphi(FloatType distance, int order)
{
   FloatType dphi = 0.0;
   FloatType sign = (distance >= 0.) ? 1.0 : -1.0;
   FloatType abs_distance = fabs(distance);

   if (order == 4) {
      FloatType abs_distance2 = abs_distance*abs_distance;
      if (abs_distance == 0.0) {
         dphi = 0.0;
      } else if (abs_distance <= 1.) {
         dphi = sign * (4.5 * abs_distance2 - 5. * abs_distance);
      } else if (abs_distance <= 2.) {
         dphi = sign * (2. - abs_distance) * (1.5 * abs_distance - 2.);
      } else {
         dphi = 0.0;
      }
   } else if (order == 6) {
      FloatType abs_distance2 = abs_distance*abs_distance;
      FloatType abs_distance3 = abs_distance2*abs_distance;
      FloatType abs_distance4 = abs_distance2*abs_distance2;
      if (abs_distance == 0.0) {
         dphi = 0.0;
      } else if (abs_distance <= 1) {
         dphi = sign * (52.*abs_distance3 - 25. * abs_distance4 + 15.*abs_distance2 - 50.*abs_distance)/12.0;
      } else if (abs_distance <= 2) {
         dphi = sign*(15*abs_distance2*abs_distance2 - 60*abs_distance2*abs_distance + 55*abs_distance2 +
               10*abs_distance4 - 96*abs_distance3 + 260*abs_distance2 - 210*abs_distance + 10)/ 24.0;
      } else if (abs_distance <= 3) {
         dphi = - sign *(abs_distance - 3)*(5*abs_distance3 - 37*abs_distance2 +
               84*abs_distance - 58)/24.0;
      } else {
         dphi = 0.0;
      }

   } else if (order == 8) {
      FloatType abs_distance2 = abs_distance*abs_distance;
      FloatType abs_distance3 = abs_distance2*abs_distance;
      FloatType abs_distance4 = abs_distance2*abs_distance2;
      FloatType abs_distance5 = abs_distance4*abs_distance;
      FloatType abs_distance6 = abs_distance5*abs_distance;
      if (abs_distance == 0.0) {
         dphi = 0.0;
      } else if (abs_distance <= 1) {
         dphi = sign*(7*abs_distance6 + 42*abs_distance4*abs_distance2 - 134*abs_distance4*abs_distance - 35*abs_distance4 -
               16*abs_distance2*abs_distance3 - 140*abs_distance2*abs_distance2 + 604*abs_distance2*abs_distance + 28*abs_distance2 +
               40*abs_distance3 + 56*abs_distance2 - 560*abs_distance)/144.0;
      } else if (abs_distance <= 2) {
         dphi = sign * (126*abs_distance4*abs_distance - 21*abs_distance4*abs_distance2 - 182*abs_distance4 -
               28*abs_distance2*abs_distance4 + 300*abs_distance2*abs_distance3 - 1001*abs_distance2*abs_distance2 +
               990*abs_distance2*abs_distance + 154*abs_distance2 + 24*abs_distance5 - 182*abs_distance4 +
               270*abs_distance3 + 602*abs_distance2 - 1260*abs_distance + 28)/240.0;
      } else if (abs_distance <= 3) {
         dphi = sign * (35*abs_distance2*abs_distance4 - 420*abs_distance2*abs_distance3 +
               1785*abs_distance2*abs_distance2 - 3150*abs_distance2*abs_distance + 1918*abs_distance2 +
               14*abs_distance6 - 330*abs_distance5 + 2660*abs_distance4 -
               9590*abs_distance3 + 15806*abs_distance2 - 9940*abs_distance + 756)/720.0;
      } else if (abs_distance <= 4) {
         dphi = -sign *(abs_distance - 4)*(7*abs_distance5 - 122*abs_distance4 +
               807*abs_distance3 - 2512*abs_distance2 + 3644*abs_distance - 1944)/720.0;
      } else {
         dphi = 0.0;
      }

   } else if (order == 10) {
      FloatType abs_distance2 = abs_distance*abs_distance;
      FloatType abs_distance3 = abs_distance2*abs_distance;
      FloatType abs_distance4 = abs_distance2*abs_distance2;
      FloatType abs_distance5 = abs_distance4*abs_distance;
      FloatType abs_distance6 = abs_distance5*abs_distance;
      FloatType abs_distance7 = abs_distance6*abs_distance;
      FloatType abs_distance8 = abs_distance7*abs_distance;
      if (abs_distance == 0.0) {
         dphi = 0.0;
      } else if (abs_distance <= 1) {
         dphi = sign*(328.*abs_distance7 - 81.*abs_distance8 +
               882.*abs_distance6 - 3644.*abs_distance5 -
               441.*abs_distance4 - 280.*abs_distance5 - 1764.*abs_distance4 + 12026.*abs_distance3 +
               324.*abs_distance2 + 490.*abs_distance3 + 648.*abs_distance2 - 10792.*abs_distance)/2880.0;
      } else if (abs_distance <= 2) {
         dphi = sign * (9.*abs_distance8 - 72.*abs_distance7 + 141.*abs_distance6 +
               18.*abs_distance8 - 236.*abs_distance7 + 963.*abs_distance6 -
               1046.*abs_distance5 - 687.*abs_distance4 - 20.*abs_distance7 + 156.*abs_distance6 +
               168.*abs_distance5 - 3522.*abs_distance4 + 6382.*abs_distance3 + 474.*abs_distance2 +
               50.*abs_distance5 - 516.*abs_distance4 + 1262.*abs_distance3 + 1596.*abs_distance2 -
               6344.*abs_distance + 72.)/1440.0;
      } else if (abs_distance <= 3) {
         dphi = sign*(720.*abs_distance7 - 45.*abs_distance8 - 4185.*abs_distance6 +
               10440.*abs_distance5 - 9396.*abs_distance4 - 36.*abs_distance8 + 870.*abs_distance7 -
               7965.*abs_distance6 + 34540.*abs_distance5 - 70389.*abs_distance4 +
               51440.*abs_distance3 + 6012.*abs_distance2 + 50.*abs_distance7 - 954.*abs_distance6 +
               6680.*abs_distance5 - 19440.*abs_distance4 + 11140.*abs_distance3 + 49014.*abs_distance2 -
               69080.*abs_distance + 3384.)/10080.0;
      } else if (abs_distance <= 4) {
         dphi = sign *(63.*abs_distance8 - 1512.*abs_distance7 + 14490.*abs_distance6 -
               70560.*abs_distance5 + 182763.*abs_distance4 - 236376.*abs_distance3 +
               117612.*abs_distance2 + 18.*abs_distance8 - 784.*abs_distance7 + 12600.*abs_distance6 -
               101556.*abs_distance5 + 451962.*abs_distance4 - 1121316.*abs_distance3 +
               1451628.*abs_distance2 - 795368.*abs_distance + 71856.)/40320.0;
      } else if (abs_distance <= 5) {
         dphi = -sign*(abs_distance - 5.)*(9.*abs_distance7 - 283.*abs_distance6 +
               3667.*abs_distance5 - 25261.*abs_distance4 + 99340.*abs_distance3 -
               221416.*abs_distance2 + 256552.*abs_distance - 117648.)/40320.0;
      } else {
         dphi = 0.0;
      }
   }
   return dphi;
}



void msm::functions::init()
{ 
   memory::create(msm::functions::gcons,7,7);
   gcons[2][0] = 15.0 / 8.0;
   gcons[2][1] = -5.0 / 4.0;
   gcons[2][2] = 3.0 / 8.0;
   gcons[3][0] = 35.0 / 16.0;
   gcons[3][1] = -35.0 / 16.0;
   gcons[3][2] = 21.0 / 16.0;
   gcons[3][3] = -5.0 / 16.0;
   gcons[4][0] = 315.0 / 128.0;
   gcons[4][1] = -105.0 / 32.0;
   gcons[4][2] = 189.0 / 64.0;
   gcons[4][3] = -45.0 / 32.0;
   gcons[4][4] = 35.0 / 128.0;
   gcons[5][0] = 693.0 / 256.0;
   gcons[5][1] = -1155.0 / 256.0;
   gcons[5][2] = 693.0 / 128.0;
   gcons[5][3] = -495.0 / 128.0;
   gcons[5][4] = 385.0 / 256.0;
   gcons[5][5] = -63.0 / 256.0;
   gcons[6][0] = 3003.0 / 1024.0;
   gcons[6][1] = -3003.0 / 512.0;
   gcons[6][2] = 9009.0 / 1024.0;
   gcons[6][3] = -2145.0 / 256.0;
   gcons[6][4] = 5005.0 / 1024.0;
   gcons[6][5] = -819.0 / 512.0;
   gcons[6][6] = 231.0 / 1024.0;

   memory::create(msm::functions::dgcons,7,6);
   dgcons[2][0] = -5.0 / 4.0 * 12.;
   dgcons[2][1] = 3.0 / 8.0 * 24.;
   dgcons[3][0] = -35.0 / 16.0 * 12.;
   dgcons[3][1] = 21.0 / 16.0 * 24.;
   dgcons[3][2] = -5.0 / 16.0 * 36.;
   dgcons[4][0] = -105.0 / 32.0 * 12.;
   dgcons[4][1] = 189.0 / 64.0 * 24.;
   dgcons[4][2] = -45.0 / 32.0 * 36.;
   dgcons[4][3] = 35.0 / 128.0 * 48.;
   dgcons[5][0] = -1155.0 / 256.0 * 12.;
   dgcons[5][1] = 693.0 / 128.0 * 24.;
   dgcons[5][2] = -495.0 / 128.0 * 36.;
   dgcons[5][3] = 385.0 / 256.0 * 48.;
   dgcons[5][4] = -63.0 / 256.0 * 60.;
   dgcons[6][0] = -3003.0 / 512.0 * 12.;
   dgcons[6][1] = 9009.0 / 1024.0 * 24.;
   dgcons[6][2] = -2145.0 / 256.0 * 36.;
   dgcons[6][3] = 5005.0 / 1024.0 * 48.;
   dgcons[6][4] = -819.0 / 512.0 * 60.;
   dgcons[6][5] = 231.0 / 1024.0 * 72.;
}
