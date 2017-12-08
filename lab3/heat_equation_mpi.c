#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <png.h>
#include <stdbool.h>
#include <mpi.h>

#define X 0                        //-
#define Y 1                        //-
#define HEQ_INIT_TAG 2
#define HEQ_FINISH_TAG       3        //-
#define CURRENT_TIMESTEP  4        //-


int           rank_global;         // The current MPI rank in the global communicator.
int           rank;                // The current MPI rank in the local communicator.
int           num_tasks;           //- The number of processes in the current local communicator
int           num_tasks_global;    //- The global number of processe

MPI_Datatype  filetype_temp;       //-
MPI_Datatype  filetype_mat;       //-
MPI_Group basegroup;               //-
MPI_Group     group;               //- A group of all workers, excluding the master process
MPI_Comm      cart_comm;           //- Communicator for the cartesian grid
MPI_Comm      workerscomm;         //- Communucator for the workers
MPI_File      file;                // A shared file pointer
MPI_Datatype  memtype_temp;             // A new created type for the inner and outer data including the ghost layer
MPI_Datatype  memtype_mat;             // A new created type for the inner and outer data including the ghost layer


#ifdef __ORDER_LITTLE_ENDIAN__
inline static double double_swap( double d )
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    union
    {
       double d;
       char b[8];
    } dat1, dat2;

    dat1.d = d;
    dat2.b[0] = dat1.b[7];
    dat2.b[1] = dat1.b[6];
    dat2.b[2] = dat1.b[5];
    dat2.b[3] = dat1.b[4];
    dat2.b[4] = dat1.b[3];
    dat2.b[5] = dat1.b[2];
    dat2.b[6] = dat1.b[1];
    dat2.b[7] = dat1.b[0];
    return dat2.d;
#else
    return d;
#endif
}
#else
#error "No endianness macro defined"
#endif


#define calcIndex(width, x,y)  ((y)*(width) + (x))
const unsigned char color_black[3] = {0,0,0};
const unsigned char color_white[3] = {255,255,255};
const unsigned char color_red[3]   = {255,0,0};
const unsigned char color_blue[3]  = {0,0,255};
const unsigned char color_green[3] = {0,255,0};
const unsigned char color_yellow[3]= {255,255,0};

#define BOUNDARY_NONE 0
#define BOUNDARY_CONST_TEMP 1
#define BOUNDARY_NEUMANN 2
#define ICE 0
#define BEER 1
#define GLASS 2

double dt = 1.0;
double dxR2 = 1.0;

typedef struct Material {
  double lambda;    // thermal conductivity W/(mK)
  double c;         // heat capacity  J/(m^3K)
  double cR;        // inverse heat capacity
  double initT;     // initial temperature Degree C
  int boundary;     // boundary type (BOUNDARY_NONE is default)
  char * name;      // string name of the material
  unsigned char color[3]; // png color to be assigned by the init_filling function
} Material;

void Material_init (Material * self, const char * name, double initT,
                    const unsigned char * color, double lambda, double c, int boundary) {
  self->name = calloc (strlen (name) + 1, sizeof(char));
  strcpy (self->name, name);
  self->color[0] = color[0];
  self->color[1] = color[1];
  self->color[2] = color[2];
  self->lambda = lambda;
  self->c = c;
  self->cR = 1 / c;
  self->initT = initT;
  self->boundary = boundary;
}

void myexit (const char * s, ...) {
  va_list args;
  va_start(args, s);
  vprintf (s, args);
  printf ("\n");
  va_end(args);
  abort ();
}

char vtk_header[2048];
void create_vtk_header (char * header, int width, int height, int timestep) {
  char buffer[1024];
  header[0] = '\0';
  strcat (header, "# vtk DataFile Version 3.0\n");
  snprintf (buffer, sizeof(buffer), "Heat equation timestep %d \n", timestep);
  strcat (header, buffer);
  strcat (header, "BINARY\n");
  strcat (header, "DATASET STRUCTURED_POINTS\n");
  snprintf (buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat (header, buffer);
  strcat (header, "SPACING 1.0 1.0 1.0\n");
  strcat (header, "ORIGIN 0 0 0\n");
  snprintf (buffer, sizeof(buffer), "POINT_DATA %ld\n", width * height);
  strcat (header, buffer);
  strcat (header, "SCALARS data double 1\n");
  strcat (header, "LOOKUP_TABLE default\n");
}

png_bytep * row_pointers;
png_byte color_type;
void read_png_file (char* file_name, int *pwidth, int *pheight) {
  int x, y;
  png_byte bit_depth;

  png_structp png_ptr;
  png_infop info_ptr;
  int number_of_passes;

  char header[8];    // 8 is the maximum size that can be checked

  /* open file and test for it being a png */
  FILE *fp = fopen (file_name, "rb");
  if (!fp)
    myexit ("[read_png_file] File %s could not be opened for reading", file_name);
  fread (header, 1, 8, fp);
  if (png_sig_cmp (header, 0, 8))
    myexit ("[read_png_file] File %s is not recognized as a PNG file", file_name);

  /* initialize stuff */
  png_ptr = png_create_read_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr)
    myexit ("[read_png_file] png_create_read_struct failed");

  info_ptr = png_create_info_struct (png_ptr);
  if (!info_ptr)
    myexit ("[read_png_file] png_create_info_struct failed");

  if (setjmp(png_jmpbuf(png_ptr)))
    myexit ("[read_png_file] Error during init_io");

  png_init_io (png_ptr, fp);
  png_set_sig_bytes (png_ptr, 8);

  png_read_info (png_ptr, info_ptr);

  (*pwidth) = png_get_image_width (png_ptr, info_ptr);
  (*pheight) = png_get_image_height (png_ptr, info_ptr);
  int width = (*pwidth);
  int height = (*pheight);

  printf ("PNG resolution (x times y) = (%d x %d)\n", width, height);

  color_type = png_get_color_type (png_ptr, info_ptr);
  bit_depth = png_get_bit_depth (png_ptr, info_ptr);
  if (bit_depth != 8) {
    myexit ("%s", "Only png files with bit depth 8 supported");
  }

  if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
    printf ("PNG color-type = %s\n", "RGB_ALPHA");
  }
  else if (color_type == PNG_COLOR_TYPE_RGB) {
    printf ("PNG color-type = %s\n", "RGB");
  }
  else {
    myexit ("%s", "Only png files with RGB and RGBA color-type supported!");
  }
  //printf("color_type = %d, bit_depth = %d\n",(int)color_type,(int) bit_depth);
  number_of_passes = png_set_interlace_handling (png_ptr);
  png_read_update_info (png_ptr, info_ptr);

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr)))
    myexit ("[read_png_file] Error during read_image");

  row_pointers = (png_bytep*) malloc (sizeof(png_bytep) * height);
  for (y = 0; y < height; y++)
    row_pointers[y] = (png_byte*) malloc (png_get_rowbytes (png_ptr, info_ptr));

  png_read_image (png_ptr, row_pointers);

  fclose (fp);
}

unsigned char * get_color_from_png_data (int png_x, int png_y, int png_height) {
  png_byte* color_row = row_pointers[png_height - png_y - 1];
  int offset = 0;
  if (color_type == PNG_COLOR_TYPE_RGB) {
    offset = 3;
  }
  else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
    offset = 4;
  }
  return &(color_row[png_x * offset]);
}


void write_field (double* currentfield, double * oldfield, int * gsizes, int * lsizes, int timestep) {

  if (timestep == 0){
     if(rank == 0) {
       mkdir("./heq/", 0777);
     }
     create_vtk_header (vtk_header, gsizes[X], gsizes[Y], timestep);
   }

   printf ("rank %d writing timestep %d\n", rank_global, timestep);
   char filename[1024];
   snprintf (filename, 1024, "./heq/heq-%05d.vtk", timestep);
   MPI_Offset header_offset = (MPI_Offset)strlen(vtk_header);

   /* TODO Create a new file handle for collective I/O
    *      Use the global 'file' variable.
    * HINT: Do not forget to swap the byte-ordering for the double field before writing to the file.
    *       Use therefore the 'double_swap' function.
    */

   /* TODO Set the file view for the file handle using collective I/O
    *
    */

   /* TODO Write the data using collective I/O
    *
    */

   /* TODO Close the file handle.
    *
    */

}

void evolve (double* current_tempfield, double* new_tempfield, int * materialfield,
             Material * mats, int width, int height) {
  double flow[4];
  int mat_ids[4];
  int mat_center;
  double temp_center;
  int i_center;
  int i_neighbors[4];
  // TODO implement calculation of average beer temperature and print out in the console
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i_center = calcIndex(width, x, y);
      i_neighbors[0] = calcIndex(width, x - 1, y);
      i_neighbors[1] = calcIndex(width, x + 1, y);
      i_neighbors[2] = calcIndex(width, x, y - 1);
      i_neighbors[3] = calcIndex(width, x, y + 1);
      mat_center = materialfield[i_center];
      if (mats[mat_center].boundary == BOUNDARY_NONE) {
        mat_ids[0] = materialfield[i_neighbors[0]];
        mat_ids[1] = materialfield[i_neighbors[1]];
        mat_ids[2] = materialfield[i_neighbors[2]];
        mat_ids[3] = materialfield[i_neighbors[3]];

        temp_center = current_tempfield[i_center];
        flow[0] = ((mats[mat_ids[0]].lambda + mats[mat_center].lambda) * 0.5) * (temp_center - current_tempfield[i_neighbors[0]]);
        flow[1] = ((mats[mat_ids[1]].lambda + mats[mat_center].lambda) * 0.5) * (current_tempfield[i_neighbors[1]] - temp_center);
        flow[2] = ((mats[mat_ids[2]].lambda + mats[mat_center].lambda) * 0.5) * (temp_center - current_tempfield[i_neighbors[2]]);
        flow[3] = ((mats[mat_ids[3]].lambda + mats[mat_center].lambda) * 0.5) * (current_tempfield[i_neighbors[3]] - temp_center);
        const double delta_temp = dt * dxR2 * mats[mat_center].cR      * (-flow[0] + flow[1] - flow[2] + flow[3]);
        new_tempfield[i_center] = temp_center + delta_temp;
      }
    }
  }
}

void init_filling (int *mat_field, double *tempfield, Material * mats, int n_mats, int width,
                   int height) {
  int i;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      i = calcIndex(width, x, y);
      int material_value = 0;
      double T = 0.0;
      double min_diff = 1e99;
      double diff;
      unsigned char * color = get_color_from_png_data (x, y, height);
      for (int ic = 0; ic < n_mats; ic++) {
        diff = fabs ((double) (color[0] - mats[ic].color[0]));
        diff += fabs ((double) (color[1] - mats[ic].color[1]));
        diff += fabs ((double) (color[2] - mats[ic].color[2]));
        if (diff < min_diff) {
          material_value = ic;
          min_diff = diff;
        }
      }
      mat_field[i] = material_value;
      tempfield[i] = mats[material_value].initT;
    }
  }
}

void heat_equation (int num_timesteps, int px, int py, int process_numX, int process_numY) {
  MPI_Status status;
  Material materials[3];
  Material_init (&materials[ICE], "ice", 0.0, color_white, 0.0, 0.0, BOUNDARY_CONST_TEMP);
  Material_init (&materials[BEER], "beer", 25.0, color_yellow, 1.0, 5.0, BOUNDARY_NONE);
  Material_init (&materials[GLASS], "glass", 18.0, color_black, 0.1, 2.0, BOUNDARY_NONE);
  /* TODO read image file with the master node and distribute parts to the worker nodes
   *
   */
}

int main (int c, char **v) {

  MPI_Init(&c, &v);

  int width = 0, height = 0, num_timesteps;
  int process_numX;
  int process_numY;

  if (c == 4) {
    num_timesteps = atoi (v[1]); ///< read timesteps

    if (width <= 0) {
      width = 32; ///< default width
    }
    if (height <= 0) {
      height = 32; ///< default height
    }
    process_numX = atoi (v[2]); ///< read number of processes in X
    process_numY = atoi (v[3]); ///< read number of processes in Y
  }
  else {
   myexit("Too less arguments");
  }

  // TODO get the global rank of the process and save it to rank_global
  // TODO get the number of processes and save it to num_tasks_global variable

  /* Abort if the number of processes does not match with the given configuration.
   */
  if (num_tasks_global != (process_numX*process_numY+1)) {
    if (rank_global == 0) {
      fprintf(stderr, "ERROR: %d MPI processes needed.\n", process_numX*process_numY+1);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }


  MPI_Comm_group(MPI_COMM_WORLD,&basegroup);
  int ranks_to_exclude[1] = { 0 };
  MPI_Group_excl(basegroup, 1, ranks_to_exclude, &group);

  // communicator of all workers
  MPI_Comm_create(MPI_COMM_WORLD, group, &workerscomm);

  if (rank_global == 0) {

    printf("[MASTER] with rank: %d is waiting...\n", rank_global);
    MPI_Status status;
    char       buffer[4096];
    //int *t= calloc(1,sizeof(int));
    int t=0;
    int        count;
    bool       running = true;
    heat_equation(num_timesteps, -1, -1, process_numX, process_numY);
    while (running) {
      // probe and wait for any MPI message
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      switch (status.MPI_TAG) {
        case 1:
          MPI_Get_count(&status, MPI_CHAR, &count);
          MPI_Recv(buffer, count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          break;
        case HEQ_FINISH_TAG:
          running = false;
          break;
        case CURRENT_TIMESTEP:
          MPI_Get_count(&status, MPI_INT, &count);
          MPI_Recv(&t, count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          printf("-----    timestep %d mpi_source %ld\n",t,status.MPI_SOURCE);
          break;
        default:
          fprintf(stderr, "[MASTER]: Unexpected MPI_TAG received (%d).\n", status.MPI_TAG);
          running = false;
      }
    }
  } else {
    /* TODO Create a new cartesian communicator of the worker communicator and get the information.
    */

    if (1 == rank) {
      long buffer = 0;
      MPI_Send(&buffer, 1, MPI_LONG, 0, HEQ_FINISH_TAG, MPI_COMM_WORLD);
    }

  }
  MPI_Finalize();
}
