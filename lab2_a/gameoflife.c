#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

//OPTIONAL: comment this out for console output
//#define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

void myexit(const char *s, ...)
{
  va_list args;
  va_start(args, s);
  vprintf(s, args);
  printf("\n");
  va_end(args);
  abort();
}

char vtk_header[2048];
void create_vtk_header(char *header, int width, int height, int timestep)
{
  char buffer[1024];
  header[0] = '\0';
  strcat(header, "# vtk DataFile Version 3.0\n");
  snprintf(buffer, sizeof(buffer), "Gameoflife timestep %d \n", timestep);
  strcat(header, buffer);
  strcat(header, "BINARY\n");
  strcat(header, "DATASET STRUCTURED_POINTS\n");
  snprintf(buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat(header, buffer);
  strcat(header, "SPACING 1.0 1.0 1.0\n");
  strcat(header, "ORIGIN 0 0 0\n");
  snprintf(buffer, sizeof(buffer), "POINT_DATA %ld\n", width * height);
  strcat(header, buffer);
  strcat(header, "SCALARS data char 1\n");
  strcat(header, "LOOKUP_TABLE default\n");
}

void write_vtk_data(FILE *f, char *data, int length)
{
  if (fwrite(data, sizeof(char), length, f) != length)
  {
    myexit("Could not write vtk-Data");
  }
}

void write_field(char *currentfield, int width, int height, int timestep)
{
#ifdef CONSOLE_OUTPUT
  printf("\033[H");
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
      printf(ALIVE == currentfield[calcIndex(width, x, y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
  printf("\ntimestep=%d", timestep);
  usleep(80000);
#else
  if (timestep == 0)
  {
    mkdir("./gol/", 0777);
    create_vtk_header(vtk_header, width, height, timestep);
  }
  printf("writing timestep %d\n", timestep);
  FILE *fp; // The current file handle.
  char filename[1024];
  snprintf(filename, 1024, "./gol/gol-%05d.vtk", timestep);
  fp = fopen(filename, "w");
  write_vtk_data(fp, vtk_header, strlen(vtk_header));
  write_vtk_data(fp, currentfield, width * height);
  fclose(fp);
  printf("finished writing timestep %d\n", timestep);
#endif
}

int countLifingsPeriodic(char *currentfield, int x, int y, int w, int h)
{
  int n = 0;
  for (int y1 = y - 1; y1 <= y + 1; y1++)
  {
    for (int x1 = x - 1; x1 <= x + 1; x1++)
    {
      if (currentfield[calcIndex(w, (x1 + w) % w, (y1 + h) % h)])
      {
        n++;
      }
    }
  }
  return n;
}

void evolve(char *currentfield, char *newfield, int width, int height)
{
  // TODO traverse through each voxel and implement game of live logic
  // HINT: avoid boundaries
#pragma omp parallel for
  for (int y = 1; y < height - 1; y++)
  {
    for (int x = 1; x < width - 1; x++)
    {
      int n = countLifingsPeriodic(currentfield, x, y, width, height);
      if (currentfield[calcIndex(width, x, y)])
      {
        n--;
      }
      newfield[calcIndex(width, x, y)] = (n == 3 || (n == 2 && currentfield[calcIndex(width, x, y)]));
    }
  }
}

void filling_random(char *currentfield, int width, int height)
{
  int i;
  for (int y = 1; y < height - 1; y++)
  {
    for (int x = 1; x < width - 1; x++)
    {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner(char *currentfield, int width, int height)
{
  currentfield[calcIndex(width, width / 2 + 0, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 1, height / 2 + 2)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 0)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 2)] = ALIVE;
}

void fill_boundaries(char *currentField, int width, int height)
{
  for (int i = 1; i < width - 1; i++)
  {
    currentField[calcIndex(width, i, 0)] = currentField[calcIndex(width, i, height - 2)];
    currentField[calcIndex(width, i, height - 1)] = currentField[calcIndex(width, i, 1)];
  }
  for (int i = 1; i < height - 1; i++)
  {
    currentField[calcIndex(width, 0, i)] = currentField[calcIndex(width, width - 2, i)];
    currentField[calcIndex(width, width - 1, i)] = currentField[calcIndex(width, 1, i)];
  }
  currentField[0] = currentField[calcIndex(width, width - 2, height - 2)];
  currentField[calcIndex(width, width - 1, height - 1)] = currentField[calcIndex(width, 1, 1)];
  currentField[calcIndex(width, width - 1, 0)] = currentField[calcIndex(width, 1, height - 2)];
  currentField[calcIndex(width, 0, height - 1)] = currentField[calcIndex(width, width - 2, 1)];
}

void game(int width, int height, int num_timesteps)
{
  char *currentfield = calloc(width * height, sizeof(char));
  char *newfield = calloc(width * height, sizeof(char));

  // TODO 1: use your favorite filling
  //filling_random(currentfield, width, height);
  filling_runner(currentfield, width, height);

  int time = 0;
  write_field(currentfield, width, height, time);
  fill_boundaries(currentfield, width, height);
  for (time = 1; time <= num_timesteps; time++)
  {
    // TODO 2: implement evolve function(see above)
    evolve(currentfield, newfield, width, height);
    //write_field(newfield, width, height, time);

    // TODO 3: implement SWAP of the fields
    char *temp = currentfield;
    currentfield = newfield;
    newfield = temp;

    // TODO 4: implement periodic boundary condition
    fill_boundaries(currentfield, width, height);
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v)
{
  int width = 0, height = 0, num_timesteps;
  if (c == 4)
  {
    width = atoi(v[1]) + 2;     ///< read width + 2 boundary cells(low x, high x)
    height = atoi(v[2]) + 2;    ///< read height + 2 boundary cells(low y, high y)
    num_timesteps = atoi(v[3]); ///< read timesteps
    if (width <= 0)
    {
      width = 32; ///< default width
    }
    if (height <= 0)
    {
      height = 32; ///< default height
    }
    game(width, height, num_timesteps);
  }
  else
  {
    myexit("Too less arguments");
  }
}
