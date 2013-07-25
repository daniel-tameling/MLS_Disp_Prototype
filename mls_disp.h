
#ifndef __MSM_H
#define __MSM_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

#include "util.h"


namespace msm{
    class Grid;
    class Config;
    class ShortRange;
    class LongRange;
    class Cell;
}

/// sets all three pointers with NULL //TODO auslagern in util
template<class T>
void reset_ptr3(T* (&ptr)[3])
{
    ptr[0] = NULL;
    ptr[1] = NULL;
    ptr[2] = NULL;
}

/// allocates all three pointers //TODO auslagern in util
template<class T>
void alloc_ptr3(T* (&ptr)[3], int num_elements)
{
    ptr[0] = new T[num_elements];
    ptr[1] = new T[num_elements];
    ptr[2] = new T[num_elements];
}

/// Frees allocated data of the given pointer
template<class T>
void free_ptr3(T* (&ptr)[3])
{
    if(ptr[0] != NULL)
        delete [] ptr[0];
    if(ptr[1] != NULL)
        delete [] ptr[1];
    if(ptr[2] != NULL)
        delete [] ptr[2];
}



namespace msm{
    /**
     * \namespace msm 
     * \brief acts like a main for a MSM simulation
     **/

    /// All data related to the configuration of a MSM simulation is encapsulated in this class
    class Config{
       public:
          /// read lammps config file 
          Config() {};
          Config(std::string input_filename, std::string config_filename);
          ~Config();

          FloatType **epsilon;
          FloatType **sigma;
          FloatType *sqrt_epsilon; //auxiliary array //TODO may not be necessary
          FloatType *sigma3; //auxiliary array  //TODO may not be necessary
          FloatType *particle_mass;
          FloatType domain_extend[3];
          FloatType dt;
          int num_particle_types;
          int output_frequency;
          int num_particles;
          int num_iterations;
          int num_cells[3];
          int num_grids;
          int num_gridpoints[3];
          int interpolation_order;
          int split_order;
          int stencil_size; /// stincil size which is used by anter-/ inter-polation and restriction and prolongation.
          FloatType grid_spacing[3];
          FloatType cutoff, cutoff_sqr;
    };


    /// Starts the simulation
    void run(msm::Config* config, ShortRange* short_range, LongRange* long_range);
        
namespace solver{
   /** 
    * \namespace msm::solver 
    * \brief Integrates Newtons equation of motion by using the velocity verlet scheme
    */

    /// integrates All given particles in time.
    void verlet(  const FloatType dt,
            const int num_particles,
            int* particle_type,
            FloatType* particle_mass,
            FloatType** particle_force,
            FloatType** particle_force_old, 
            FloatType** particle_vel, 
            FloatType** particle_pos);
    
    /// unit test for verlet()
    void test_verlet();
}

namespace functions{
    /**
     * \namespace msm::functions
     * \brief provides the necessary mathematical functions used by msm
     **/

   /// gamma function used to bring the potential to zero as the distance approaches the cutoff radius.
   FloatType gamma( FloatType distance_sqr, int split_order);
   /// first derivative of gamma devided by distance
   FloatType dgamma( FloatType distance_sqr, int split_order);

   /// computes the phi function used for the grid approximation using an interpolation order of 2
   FloatType phi(FloatType distance, int order);
   /// computes the derivative of phi used for the grid approximation using an interpolation order of 2
   FloatType dphi(FloatType distance, int order);

   void init();
   extern FloatType** gcons;
   extern FloatType** dgcons;

}

/**
 * \class LongRange
 * \brief Stores all data and provides all functions required for the long-range part of a msm simulation.
 */
class LongRange{
    
    public:
        /// initializes the longRange part
        LongRange(msm::Config* config);

        /// frees alocated memory (i.e. grids). This is the counterpart to init.
        ~LongRange();

        /// computes one step of the long range computation and stores the forces to the Cells
        void computeForces(int num_cells[3], Cell*** cells); 

        /// get potential Energy. WARNING: it does not recompute it.
        FloatType getPotentialEnergy() { return this->potential_energy; }

    /*-----------------------------------------------------------------------------
     *  longRange variables
     *-----------------------------------------------------------------------------*/
    private:
        /// number of grids
        int        num_grids;  
        /// stores num_grids many grids. grids[0] is a FinestGrid, while girds[num_grids-1] is a CoarsestGrid.
        Grid       **grids;     
        /// potential energy
        FloatType potential_energy;
        /// precomputed phi for restriction and prolongation
        FloatType *phi_grid;
        /// indices corresponding to phi_grid
        int *phi_grid_index;
        /// precomputed gamma stencil for grid; used for calculation the potential on intermediate grids; 
        FloatType *gamma_grid;
        /// indices corresponding to gamma_grid
        int *gamma_grid_index[3];
        /// size of array gamma_grid, gamma_grid_index[0],[1] and [2]
        int gamma_grid_size;

        msm::Config* config;
};
 
/**
 * \class Grid.
 * \brief Encapsulates all data and functions associated with one particular grid level.
 **/
class Grid{
    public:

        /// Grid constructor sets all member variables to zero.
        Grid();

        /// Grid Destructor frees all allocated data
        virtual ~Grid();

        /*-----------------------------------------------------------------------------
         *  Public Functions
         *-----------------------------------------------------------------------------*/
        /// initializes all data required by the grid
        virtual void init(msm::Config* config, int level);

        ///update the potential of this grid
        virtual void direct(FloatType cutoff, FloatType* gamma_grid=0, int** gamma_grid_index=0, int gamma_grid_size=0);   

        /// sets all xi's and potentials at every gridpoint to zero
        void resetGrid();

        /// updates the xi of the dstGrid via restriction of this grid. dstGrid is one level higher than current grid (dstGrid is coarser).
        virtual void restriction(Grid* dstGrid, const int stencil_size, FloatType *phi_grid, int* phi_grid_index);

        /// updates the potential of the dstGrid vie prolongation of this grid. dstGrid is one level lower than current grid (dstGrid is finer).
        virtual void prolongation(Grid* dstGrid, const int stencil_size, FloatType* phi_grid, int* phi_grid_index);

        /*-----------------------------------------------------------------------------
         *  Getter and setter 
         *-----------------------------------------------------------------------------*/
        /// return a 3D array of gridpoint Xi 
        FloatType*** getGridpointXi() { return this->gridpoint_xi; }
        /// return a 3D array of gridpoint potential
        FloatType*** getGridpointPotential() { return this->gridpoint_potential; }
        /// return a 3D array of gridpoint potential
        FloatType*** getGridpointForceX() { return this->gridpoint_force_x; }
        /// return a 3D array of gridpoint potential
        FloatType*** getGridpointForceY() { return this->gridpoint_force_y; }
        /// return a 3D array of gridpoint potential
        FloatType*** getGridpointForceZ() { return this->gridpoint_force_z; }

        /// return an array of size 3 containing the number of gridpoints in each dimension.
        int* getNumGridpoints() { return this->num_gridpoints; }

    protected:

        /*-----------------------------------------------------------------------------
         *  Private Variables
         *-----------------------------------------------------------------------------*/
        FloatType    spacing_grid[3]; /// spacing in each dimension of the finest grid. [0]=x-dim; [1]=y-dim; [2]=z-dim; Intermediate grids only use these values within init().
        int          num_gridpoints[3];      /// number of grid points in each dimension. [0]=x-dim; [1]=y-dim; [2]=z-dim; 
        int          level;                  /// gird level (e.g. 0 corresponds to the finest grid, 1 is coarser,...)
        FloatType    direct_aux;                  /// auxiliary variable for direct part
        FloatType*** gridpoint_xi;       /// particles on grid  of size dim[0] * dim[1] * dim[2]
        FloatType*** gridpoint_potential;    /// potential on grid of size dim[0] * dim[1] * dim[2]
        FloatType*** gridpoint_force_x;       /// x-component of force on grid of size dim[0] * dim[1] * dim[2]
        FloatType*** gridpoint_force_y;       /// y-component of force on grid of size dim[0] * dim[1] * dim[2]
        FloatType*** gridpoint_force_z;       /// z-component of force on grid of size dim[0] * dim[1] * dim[2]
        int          interpolation_order;    /// this value determins the stencil_size = 2 * interpolation_order + 1
        FloatType*** direct_stencil;         /// precomputed stencil used by the direct() function of the coarsest grid
        int          direct_periodic_int_mult_of_num_gridpoints[3]; /// auxiliary variable used in direct() for the periodic boundary computation in each dimension

};

/**
 * \class FinestGrid
 * \brief This class represents the finest grid (i.e. lowest level). It is resposible for interacting with the particles.
 **/
class FinestGrid : public Grid{
    public:
        /// TODO: implement restriction, direct, prolongation.
         
        /// find lowest left front gridpoint that is affected by the given particle
        void find_lowest_gridpoint_idx(FloatType** particle_pos, const int particle_idx, int (&lower_grid_idx)[3]);

        /// map all xi of the particles to gridpoints
        void anterpolate(Cell*** cells, const int num_cells[3], msm::Config* config);

        /// update all xi of all particle by accounting for the long-range part computed on the grid
        FloatType interpolate(Cell*** cells, const int num_cells[3], msm::Config* config);

};

/**
 * \class CoarsestGrid
 * \brief This Grid does not propagate its potential to a higher level Grid but computes an all-to-all gridPoint computation.
 **/
class CoarsestGrid : public Grid{

    public:
        virtual void init(msm::Config *config, int level);
        virtual void direct(FloatType cutoff, FloatType* gamma_grid = 0, int** gamma_grid_index = 0, int gamma_grid_size = 0);
        void restriction(Grid* dstGrid, const int stencil_size, FloatType *phi_grid, int* phi_grid_index) { std::cout<< "Warning: This function must not be called by the coarest grid.\n"; exit(-1); }

    private:
        FloatType seriescalc(FloatType dist_x,FloatType dist_y,FloatType dist_z, FloatType* domain_extend, int max_iterations, FloatType cutoff,int split_order);
};       

/**
 * \class ShortRange
 * \brief Encapsulates all data and functions necessary for the computation of the short-range part of MSM. 
 **/
class ShortRange{
    public:
        ShortRange() {};
        /// initializes the ShortRange namespace.
        ShortRange( std::string  input_file_name,
                    msm::Config* config);

        /// Free allocated data
        ~ShortRange();

        ///print particle related data to file
        void print();

        /// Read particles from input file
        void read_dump(std::string input_filename, msm::Config* config);

        /// test this class for various errors (only for debugging purposes)
        void test();

        /// this method maps a given particle to the responsible cell.        
        void mapParticleToCell( 
              int id, 
              int type, 
              FloatType force[3], 
              FloatType force_old[3],
              FloatType vel[3],
              FloatType pos[3]);

        /// TODO: check for input parameters
        /// integrate system one timestep without calculating the energy
        void timeIntegrate();
        
        /// integrate system one timestep with calculation of energy which is returned
        FloatType timeIntegrateWithEnergy();
        
        /// get potential Energy. WARNING: it does not recompute it.
        FloatType getPotentialEnergy() { return this->potential_energy; }

        /// get self interaction Energy. WARNING: it does not recompute it.
        FloatType getSelfEnergy() { return this->self_energy; }

        /// calculate the kinetic energy 
        FloatType kineticEnergy();

        /// calculate the magnitude of a force between two particles
        void computeForceMag();

        /// calculate the magnitude of a force between two particles due to Lennard-Jones potential
        void computeForceMag_LJ();

        /// calculate the magnitude of a force between two particles due to the gamma term
        void computeForceMag_MSM();

        /// calculate the potential energy between two particles
        void computePot();

        /// calculate the potential energy between two particles due to Lennard-Jones potential
        void computePot_LJ();

        /// calculate the magnitude of a force between two particles due to the gamma term
        void computePot_MSM();

        /// Unit test for computeForces()
        void computeForcesUnit();

        /// computes forces for each particle which purely result from short-range interactions
        void computeForces(int symmetric = 0);

        /// computes forces and energy for each particle which purely result from short-range interactions
        FloatType computeForcesAndEnergy();

        /// get the mixed sigma and epsilon value for the particle pairing
        void computeLJparameters();

        /// integrates all particles of all cells in time
        void integrate(FloatType dt);
        
        /// moves position of a particle to the next time step
        void updatePosition(FloatType dt);

        /// moves velocities of a particle to the next time step
        void updateVelocity();

        /// when the particle has left the domain move it back with periodic boundary conditions
        /// without moving it to the new cell this is done later
        void positionBackToDomain();

        /// check whether particle moved outside of its cell; return value is a bool whether it is the case; parameter is where the new cell indices are stored
        int checkParticleOutside(int newCellInd[3]);

        /// send particle to receive buffer of different cell
        void sendToReceiveBuffer(int currentParticle, int newCellInd[3]);

        /// integrate received particles into new cells and delete the ones that left; do this in a way that the array remains in a compressed form
        void updateCells();
        
        /// maps each particle to the domain. 
        void periodicBoundary();

        /// maps each particle to the proper cell
        void updateGeometricHashing();

        /// swaps the force and force_old ptr of each cell and initializes the current force with zero.
        void swapForcePtr();     

        /// combines Force contribution from longrang and shortRange
        void combineForces();

        /// returns private array num_cells
        int* getNumCells() { return this->num_cells; }

        /// returns private array of cells
        Cell*** getCells() { return this->cells; }

        /*-----------------------------------------------------------------------------
         *  Private variables
         *-----------------------------------------------------------------------------*/
    private:

        ///computes the cell-cell interaction for interior cells
        void computeForcesInteriorNonSym();
        ///computes the cell-cell interaction for interior cells but symetrically (i.e. f_ij = - f_ji)
        void computeForcesInteriorSym();
        ///computes the cell-cell interaction for boundary cells
        void computeForcesBoundary();
        ///computes the self interaction energy
        void computeSelfEnergy(msm::Config* config);

        msm::Config*        config; ///config file which stores all information related to the simulation (e.g. dt)
        Cell***            cells;            /// 3D Pointer to cells
        int                num_cells[3];      /// number of cells in each dimension. num_cells[0]:x-dim, num_cells[1]:y-dim, num_cells[2]:z-dim
        int                num_nbr_cells[3];   /// number of neighboring cells in each direction. E.g. numNbrCell[0] = 2 means that each  cell has to visit 2 cells to its left and 2 cells to its right.
        FloatType          cell_extend[3];    /// domain[i]/dim[i]
        FloatType          potential_energy; /// stores the potential energy
        FloatType          self_energy; /// stores the self interaction energy that is constant during the simulation

        int* tmp_particle_id; /// auxiliary array of all particle IDs used by updateGeometricHashing() TODO
        int* tmp_particle_type; /// auxiliary array of all particle types used by updateGeometricHashing() TODO
        FloatType** tmp_particle_pos; /// auxiliary array of all particle pos used by updateGeometricHashing() TODO
        FloatType** tmp_particle_force;  /// auxiliary array of all particle force used by updateGeometricHashing() TODO
        FloatType** tmp_particle_force_old; /// auxiliary array of all particle force old used by updateGeometricHashing() TODO
        FloatType** tmp_particle_vel; /// auxiliary array of all particle vel used by updateGeometricHashing() TODO
};


/** 
* \class Cell
* \brief A Cell stores particles which reside in the same portion of the computational domain. It's part of a domain decomposition.
*/
class Cell{
   public:


        /// Cell Constructor
        Cell();

        /// Cell Destructor
        ~Cell();

        /// test the cell for various errors (e.g. particle has left the domain?) (only used for debuggin purposes
        void test(FloatType* particle_mass, FloatType domain_extend[3]);

        /*-----------------------------------------------------------------------------
         *  Getter and Setter
         *-----------------------------------------------------------------------------*/
        /// returns a pointer to three pointers which point to the particle position (see particle_pos).
        FloatType** getParticlePos(){ return this->particle_pos; }        
        /// returns a pointer to three pointers which point to the particle velocity (see particle_pos).
        FloatType** getParticleVel(){ return this->particle_vel; }        
        /// returns a pointer to three pointers which point to the particle force of the current timestep(see particle_pos).
        FloatType** getParticleForceCell(){ return this->particle_force_from_cells; }        
        /// returns a pointer to three pointers which point to the particle force for the grid-contribution of the current timestep(see particle_pos).
        FloatType** getParticleForceGrid(){ return this->particle_force_from_grids; }        
        /// returns a pointer to three pointers which point to the particle force of the previous timestep(see particle_pos).
        FloatType** getParticleForceOld(){ return this->particle_force_old; }        
        /// returns a pointer to an array containing the type of each particle
        int*        getParticleType(){ return this->particle_type; }        
        /// returns a pointer to an array containing the id of each particle
        int*        getParticleID(){ return this->particle_id; }        
        /// return the number of particles in this cell
        int getNumParticles(){ return this->num_particles; }        
        /// set particle count to zero (required before mapping at a new timestep)
        void resetParticleCount() { this->num_particles = 0; }

        /*-----------------------------------------------------------------------------
         *  Public Functions
         *-----------------------------------------------------------------------------*/
        /// print all particles IDs and forces. Used for DEBUGING ONLY.
        void print(std::ofstream *output);
        /// computes the self interaction energy (without the constant in front) of all particles inside the cell
        FloatType computeSelfEnergy(msm::Config* config);
        /// allocates all memory for max_particles
        void init(int max_particles);
        /// add particle to cell
        void addParticle(  
              int id,
              int type,
              FloatType force[3], 
              FloatType force_old[3],
              FloatType vel[3],
              FloatType pos[3]);
        /// computes the interaction of all particles of this cell with all particles of nbrCell.
        FloatType computeInteractionWithCellNonSym(Cell& nbr_cell, const FloatType cutoff, FloatType** epsilon, FloatType** sigma, int split_order);     
        /// computes the interaction of all particles of this cell with all particles of nbrCell for boundary cells.
        FloatType computeInteractionWithCellNonSymPeriodic(Cell& nbr_cell, const FloatType cutoff, FloatType** epsilon, FloatType** sigma, 
              int shift_coords[3],
              FloatType domain_extend[3], int split_order);     
        /// same as computeInteractionWithCellNonSym but symetric (i.e. f_ij = - f_ji). WARNING: do not interact with the cell itself (use computeInteractionWithItselfSym())
        FloatType computeInteractionWithCellSym(Cell& nbr_cell, const FloatType cutoff_sqr, FloatType** epsilon, FloatType** sigma, int split_order);
        /// computes the interaction of all particles within this cell with its self in a symmetric fashion
        FloatType computeInteractionWithItselfSym(const FloatType cutoff_sqr, FloatType** epsilon, FloatType** sigma, int split_order);
        /// sets the array particle_force to zero
        void resetForces();
        ///  swap force and force_old ptr;
        void swapForcePtr();
        /// combines Force contribution from longrang and shortRange
        void combineForces();
        /// compute kinetic energy of all particles within this cell
        FloatType kineticEnergy(FloatType* particle_mass);

   private:
        /// unit test for kineticEnergy()
        void test_kineticEnergy(FloatType* particle_mass, int num_particle_types );
        
        /*-----------------------------------------------------------------------------
         *  Private variables
         *-----------------------------------------------------------------------------*/
        int          num_particles;        /// number of particles
        int          num_particle_types;     /// total number of different particle types 
        int          max_particles;        /// maximum number of allowed particles per cell

        FloatType**   particle_pos;         /// Pointer to particle positions in Structure of Array (SoA) format. 
                                           ///       particle_pos[p][0] points to the x-coordinate of particle p 
                                           ///  particle_pos[p][1] points to the y-coordinate of particle p 
                                           ///       particle_pos[p][2] points to the z-coordinate of particle p
        FloatType**   particle_vel;         /// Particle velocities stored SoA format (see particle_pos) 
        FloatType**   particle_force_from_cells; /// Particle forces contribution of the ShortRange class. stored in SoA format (see particle_pos)
                                                   /// This array holds the final forces after a call to ShortRange::combineForces()
        FloatType**   particle_force_from_grids; /// Particle forces stored in SoA format (see particle_pos)
        FloatType**   particle_force_old;   /// Particle old forces used for Velocity Verlet (see particle_pos)
        int*         particle_type;           /// Particle type, defines sigma and epsilon of LJ
        int*         particle_id;           /// Particle id
        int          initialized;
}; //end class Cell


}//msm

#endif
