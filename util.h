
#define sqr(x) ((x)*(x))

typedef double FloatType;

/**
 * \brief Allocates a 3D array of type <T>
 * \param   ptr   pointer to the 3D array. Allocated array has the form: ptr[num_z][num_y][num_x]. 
 * \param   num_x number of elements in x-dimension
 * \param   num_y number of elements in y-dimension
 * \param   num_z number of elements in z-dimension
 **/
template<class T>
void allocate3DArray (T*** &ptr, int num_x, int num_y, int num_z)
{
    ptr = new T** [num_z]; 
    for(int i=0;i < num_z ;++i){
        ptr[i] = new T* [num_y];
        for(int j=0;j < num_y ;++j){
            ptr[i][j] = new T[num_x];
        }
    }
}
/**
 * \brief Frees an allocated 3D array of type <T>
 * \param   ptr   pointer to the 3D array. Allocated array has the form: ptr[num_z][num_y][]. 
 * \param   num_y number of elements in y-dimension
 * \param   num_z number of elements in z-dimension
 **/
template<class T>
void free3DArray (T*** &ptr, int num_y, int num_z)
{
    for(int i=0;i < num_z ;++i){
        for(int j=0;j < num_y ;++j){
            delete [] ptr[i][j];
        }
        delete [] ptr[i];
    }
    delete [] ptr;
}

/* ----------------------------------------------------------------------
   create a 2d array
------------------------------------------------------------------------- */
namespace memory{
   //THIS FUNCTION IS TAKEN FROM memory.h
  template <typename TYPE>
    TYPE **create(TYPE **&array, int n1, int n2)
    {
      long nbytes = ((long) sizeof(TYPE)) * n1*n2;
      TYPE *data = (TYPE *) malloc(nbytes); // different from LAMMPS
      nbytes = ((long) sizeof(TYPE *)) * n1;
      array = (TYPE **) malloc(nbytes);

      long n = 0;
      for (int i = 0; i < n1; i++) {
        array[i] = &data[n];
        n += n2;
      }
      return array;
    }
}

