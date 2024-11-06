#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h> 


#define PI 3.141592654       // Constant for the value of pi
#define DECAY_FACTOR 5.0    // Factor for the distance decay in connection probability
#define INT_DIST_FACTOR 2.0  // Factor for the interdistance decay in connection probability
#define SCALE_FACTOR 1.0     // Global factor to scale the probability distribution
#define CUTOFF_RADIUS 25.0   // Maximum radius within which to consider nearby neurons
#define PACKING_FRACTION 0.5 // Maximum packing fraction of the system (area occupied / total area)
#define BRANCHES 3           // Number of branches per initial neuron


// Define a structure to represent a neuron
typedef struct
{
   double x;       // x-coordinate of the neuron
   double y;       // y-coordinate of the neuron
   double radius;  // Radius of the neuron
   int synapses;   // Number of synapses (connections)
} Neuron;


// Function prototypes
void packingFraction(int N_neurons, double L_x, double L_y, double radius);
void placeRandom(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons);
void placeLattice(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons);
void initializeNeurons(Neuron *neurons, int N_neurons, double L_x, double L_y, double radius, char mode);
void writeNeuronData(FILE *file, Neuron *neurons, int N_neurons);
double Distance(Neuron p1, Neuron p2);
bool isOverlapping(Neuron *neurons, int N_neurons, Neuron newNeuron);
float setInside(float max, float min);
bool isBetween(Neuron p1, Neuron p2, Neuron p3);
double connectionProbability(double dist, double intDistFactor);
double interDistance(Neuron p1, Neuron p2, Neuron p3);
void findNearbyParticles(Neuron *neurons, int N_neurons, int ***nearbyNeurons, int **nearbyCounts);
void freeNearbyParticles(int **nearbyNeurons, int N_neurons, int *nearbyCounts);
double*** allocate3DArray(int x, int y, int z);
void free3DArray(double ***array, int x, int y);
double** allocate2DArray(int rows, int cols);
void free2DArray(double** array, int rows);
int** allocate2DIntArray(int rows, int cols);
void free2DIntArray(int** array, int rows);


int main()
{
   int N_neurons, timeSteps, N_synapses; // Number of neurons, timeSteps for connections and synapses
   double L_x, L_y, radius;              // Physical dimensions
   char mode;                            // Mode (random|lattice)


   // Seed the random number generator
   srand(time(NULL));


   // Read the number of neurons, the dimensions of the box and neuron radii
   printf("Random or Lattice mode (R/L): ");                scanf(" %c", &mode);
   printf("Enter the number of neurons: ");                 scanf("%d", &N_neurons);
   printf("Enter the width (L_x) of the box: ");            scanf("%lf", &L_x);
   printf("Enter the height (L_y) of the box: ");           scanf("%lf", &L_y);
   printf("Enter the radius of the neurons: ");             scanf("%lf", &radius);
   printf("Enter the time steps: ");                        scanf("%d", &timeSteps);
   printf("Enter the max number of synapses per neuron: "); scanf("%d", &N_synapses);


   // Check the packing fraction to ensure it isn't too dense
   packingFraction(N_neurons, L_x, L_y, radius);


   // Allocate memory for the neurons
   Neuron *neurons = (Neuron *)calloc(N_neurons, sizeof(Neuron));
   if(!neurons)
   {
      fprintf(stderr, "Memory allocation failed\n");
      return EXIT_FAILURE;
   }


   // Initialize neuron positions
   initializeNeurons(neurons, N_neurons, L_x, L_y, radius, mode);


   // Write the neuron positions and radii to a file
   FILE *neuronFilePtr = fopen("neurons_dat", "w");
   writeNeuronData(neuronFilePtr, neurons, N_neurons);
   fclose(neuronFilePtr);


/*
+=========================================================================================================+
||                                     INITIALIZE NEURON CONNECTIONS                                     ||
+=========================================================================================================+
*/
   // Create an array of arrays (list of nearby neurons)
   int **nearbyNeurons = NULL;
   int *nearbyCounts = NULL;
   findNearbyParticles(neurons, N_neurons, &nearbyNeurons, &nearbyCounts);

   /* Create and allocate memory for a 3D array to store data for each branch
      Data: 0: N_Connections. 1: Total Distance. 2: Connected Fraction. 4: EndToEndDistance
      Array gets stored as systemData[Data][TimeSteps][Branch] */
   double ***eachBranchData = allocate3DArray(4, timeSteps, BRANCHES);
   int **neuronsList = allocate2DIntArray(BRANCHES, timeSteps);

   // Allocate arrays to keep track of each branch's progress
   int *currentNeurons = (int *)calloc(BRANCHES, sizeof(int));
   int *N_connections = (int *)calloc(BRANCHES, sizeof(int));
   double *totalDistances = (double *)calloc(BRANCHES, sizeof(double));
   double *endToEndDistances = (double *)calloc(BRANCHES, sizeof(double));
   double *connectedFractions = (double *)calloc(BRANCHES, sizeof(double));

   if(!currentNeurons || !N_connections || !totalDistances || !endToEndDistances || !connectedFractions)
   {
      fprintf(stderr, "Memory allocation failed for branch data\n");
      free(neurons);
      freeNearbyParticles(nearbyNeurons, N_neurons, nearbyCounts);
      return EXIT_FAILURE;
   }


   // Initialize all branches with the same starting random neuron
   int initialNeuron = rand() % (N_neurons + 1);
   for(int branch = 0; branch < BRANCHES; ++branch)
   {
      currentNeurons[branch] = initialNeuron;
   }
   neurons[initialNeuron].synapses = N_synapses - BRANCHES;

   // Start looking for candidates. Every branch in the same time step.
   for(int time = 0; time < timeSteps; ++time)
   {
      for(int branch = 0; branch < BRANCHES; branch++)
      {
         int currentNeuron = currentNeurons[branch];
         int nextNeuron = -1;

         // Attempt to find the next valid neuron from nearby neurons list
         for(int attempt = 0; attempt < nearbyCounts[currentNeuron]; ++attempt)
         {
            int candidateNeuron = nearbyNeurons[currentNeuron][rand() % (nearbyCounts[currentNeuron] + 1)];

            // Ensure the candidate is valid for connection
            if(candidateNeuron != initialNeuron
               && neurons[candidateNeuron].synapses < N_synapses - 1
               && neurons[currentNeuron].synapses < N_synapses)
            {
               double intDistFactor = 1.0; // Reduction factor by intermediate neurons

               // Loop through common nearby neurons for interdistance calculation
               for(int j = 0; j < nearbyCounts[currentNeuron]; j++)
               {
                  int interNeuron = nearbyNeurons[currentNeuron][j];

                  // Check if the neuron is also nearby to the candidate
                  for(int k = 0; k < nearbyCounts[candidateNeuron]; k++)
                  {
                     if(interNeuron == nearbyNeurons[candidateNeuron][k]
                        && interNeuron != currentNeuron)
                     {
                        // Calculate the interdistance
                        double interdistance = interDistance(neurons[currentNeuron], neurons[candidateNeuron], neurons[interNeuron]);

                        // Adjust the connection probability based on intermediate neurons
                        if(interdistance < radius
                           && isBetween(neurons[currentNeuron], neurons[interNeuron], neurons[candidateNeuron]))
                        {
                           if(interdistance <= (0.5 * radius)
                              && interNeuron != initialNeuron
                              && neurons[interNeuron].synapses < N_synapses - 1)
                           {
                              // The intermediate neuron becomes the new candidate
                              candidateNeuron = interNeuron;
                              intDistFactor = 1.0;
                              j = k = 0;
                           }
                           else
                           {
                              intDistFactor += interdistance / radius;
                           }
                        }
                     }
                  }
               }

               // Probability check for connection based on distance
               double dist = Distance(neurons[currentNeuron], neurons[candidateNeuron]);
               if(connectionProbability(dist, intDistFactor) > (float)rand() / RAND_MAX)
               {
                  nextNeuron = candidateNeuron;
                  break;
               }
            }
         }

         // If a connection was made, update the branch data
         if(nextNeuron != -1)
         {
            N_connections[branch]++;
            totalDistances[branch] += Distance(neurons[currentNeuron], neurons[nextNeuron]);
            connectedFractions[branch] = (double)N_connections[branch] / N_neurons;
            neurons[nextNeuron].synapses += 1;
            neurons[currentNeuron].synapses += 1;

            // Update the current neuron for the next iteration
            fprintf(stdout, "Connection made from: %d\n", currentNeuron);
            currentNeurons[branch] = nextNeuron;
         }
         else
         {
            fprintf(stderr, "No valid connection found from neuron %d\n", currentNeuron);
         }

            // Store data in multiple branches on each time step
            eachBranchData[0][time][branch] = N_connections[branch];
            eachBranchData[1][time][branch] = totalDistances[branch];
            eachBranchData[2][time][branch] = connectedFractions[branch];
            eachBranchData[3][time][branch] = Distance(neurons[initialNeuron], neurons[currentNeurons[branch]]);
            neuronsList[branch][time] = currentNeuron;

            printf("Todo guardado en el tiempo: %d  en la rama: %d\n", time, branch);
      }
   }  


   // Open the file for writing neuron connections
   FILE *connectionFilePtr = fopen("connections_dat", "w");
   if (!connectionFilePtr)
   {
      fprintf(stderr, "\nFailed to open connection data file for writing\n");
      
      // Free allocated memory
      free(neurons);
      freeNearbyParticles(nearbyNeurons, N_neurons, nearbyCounts);
      free3DArray(eachBranchData, 4, timeSteps);
      free2DIntArray(neuronsList, BRANCHES);
      free(currentNeurons);
      free(N_connections);
      free(totalDistances);
      free(endToEndDistances);
      free(connectedFractions);
   
      return EXIT_FAILURE;
   }

   // Connection format string for file output
   const char *formatConnections = "%8.4lf %8.4lf %6.2lf %10.4lf %10.0lf %15.4lf %15.4lf %15d\n";
   for(int branch = 0; branch < BRANCHES; branch++)
   {
      fprintf(connectionFilePtr, "Xlabel   Ylabel   Radius   EndToEnd   N_connections   TotalDist   ConnectedFraction   TimeStep\n");
      fprintf(connectionFilePtr, formatConnections, 
              neurons[initialNeuron].x, neurons[initialNeuron].y, neurons[initialNeuron].radius, 
              0.000000, 0.000000, 0.000000, 0.000000, 0);

      for(int time = 0; time < timeSteps; time++)
      {
         // Write connected neuron to the file
         fprintf(connectionFilePtr, formatConnections,
                 neurons[neuronsList[branch][time]].x, neurons[neuronsList[branch][time]].y, neurons[neuronsList[branch][time]].radius,
                 eachBranchData[3][time][branch], eachBranchData[0][time][branch], eachBranchData[1][time][branch],
                 eachBranchData[2][time][branch], time + 1);
      }
      //Line break to separate data
      fprintf(connectionFilePtr, "\n\n\n");
   }
   // Close connection file
   fclose(connectionFilePtr);


   // Write in a file the results for multiple branches in a single time step
   FILE *resultsFilePtr = fopen("results_dat", "w");
   if(!resultsFilePtr)
   {
      fprintf(stderr, "\nFailed to open results data file for writing\n");

      // Free allocated memory
      free(neurons);
      freeNearbyParticles(nearbyNeurons, N_neurons, nearbyCounts);
      free3DArray(eachBranchData, 4, timeSteps);
      free2DIntArray(neuronsList, BRANCHES);
      free(currentNeurons);
      free(N_connections);
      free(totalDistances);
      free(endToEndDistances);
      free(connectedFractions);


      return EXIT_FAILURE;
   }

   const char *formatResults = "%5.0lf %18.4lf %15.4lf %15d\n";
   fprintf(resultsFilePtr, "N_Connections   TotalDist   ConnectedFraction   TimeStep\n");
   fprintf(resultsFilePtr, formatResults, 0, 0.000000, 0.000000, 0);

   // Create and allocate memory for a 2D array that merges the multiple branches data into the same time step
   double** mergedBranchesData = allocate2DArray(3, timeSteps);

   // Merge the data of the 3D array into the same time step
   for(int time = 0; time < timeSteps; time++)
   {
      for(int branch = 0; branch < BRANCHES; branch++)
      {
         mergedBranchesData[0][time] += eachBranchData[0][time][branch];
         mergedBranchesData[1][time] += eachBranchData[1][time][branch];
         mergedBranchesData[2][time] += eachBranchData[2][time][branch];
      }
   }

   // Print data per time step
   for(int time = 0; time < timeSteps; time++)
   {
      fprintf(resultsFilePtr, formatResults, mergedBranchesData[0][time], mergedBranchesData[1][time], mergedBranchesData[2][time], time + 1);
   }

   fclose(resultsFilePtr);

   // Free allocated memory
   free(neurons);
   freeNearbyParticles(nearbyNeurons, N_neurons, nearbyCounts);
   // free3DArray(eachBranchData, 4, timeSteps);
   free2DArray(mergedBranchesData, 3);
   free2DIntArray(neuronsList, BRANCHES);
   free(currentNeurons);
   free(N_connections);
   free(totalDistances);
   free(endToEndDistances);
   free(connectedFractions);



   printf("Neuron data saved to neurons_dat, connection data saved to connections_dat and results_dat\n");
   return EXIT_SUCCESS;
}


/*
+=========================================================================================================+
||                                               FUNCTIONS                                               ||
+=========================================================================================================+
*/

/**************************************************************************************************************
* @brief Computes the packing fraction based on the number of neurons, box dimensions and neuron radius.
*
* @param N_neurons The total number of neurons.
* @param L_x       The width of the box.
* @param L_y       The height of the box.
* @param radius    The radius of each neuron.
*
* @return If the packing fraction is too big, the program ends.
***************************************************************************************************************/
void packingFraction(int N_neurons, double L_x, double L_y, double radius)
{
   double packingfraction = (N_neurons * PI * pow(radius, 2)) / (L_x * L_y);
   
   // Check if the packing fraction exceeds the predefined limit
   if(packingfraction >= PACKING_FRACTION)
   {
      fprintf(stderr, "Packing fraction: %.5f is too dense.\n", packingfraction);
      exit(EXIT_FAILURE);
   }
   fprintf(stdout, "Packing fraction: %.5f\n", packingfraction);
}


/**************************************************************************************************************
* @brief Places neurons inside a box in random positions avoiding overlaping.
*
* @param N_neurons The total number of neurons.
* @param L_x       The width of the box.
* @param L_y       The height of the box.
* @param radius    The radius of each neuron.
* @param neurons   Array of structures of neurons.
**************************************************************************************************************/
void placeRandom(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons)
{
   int placedParticles = 0;
      while(placedParticles < N_neurons)
      {
         Neuron newNeuron;
         newNeuron.x = setInside(L_x,radius);
         newNeuron.y = setInside(L_y,radius);
         newNeuron.radius = radius;
         newNeuron.synapses = 0;

         //Check for overlap with existing neurons
         if(!isOverlapping(neurons, placedParticles, newNeuron))
         {
            neurons[placedParticles] = newNeuron;
            placedParticles++;
         }
      }
}


/**************************************************************************************************************
* @brief Places neurons inside a box in random positions avoiding overlaping.
*
* @param N_neurons The total number of neurons.
* @param L_x       The width of the box.
* @param L_y       The height of the box.
* @param radius    The radius of each neuron.
* @param neurons   Array of structures of neurons.
**************************************************************************************************************/
void placeLattice(int N_neurons, double L_x, double L_y, double radius, Neuron *neurons)
{   
   int numRows = (int)sqrt(N_neurons);
   int numCols = (N_neurons + numRows - 1) / numRows; // Handles non-perfect squares

   double dx = (L_x - 2 * radius) / (numCols - 1);    // Spacing in the x direction
   double dy = (L_y - 2 * radius) / (numRows - 1);    // Spacing in the y direction

   double distance = sqrt(dx * dx + dy* dy);          // Handles overlapping
   if(distance < 2 * radius)
   {
      fprintf(stderr, "Error: unavoidable overlapping");
      exit(EXIT_FAILURE);
   }

   int placedParticles = 0;
   for(int i = 0; i < numRows; ++i)
   {
       for(int j = 0; j < numCols; ++j)
       {
         if(placedParticles < N_neurons)
         {
            neurons[placedParticles].x = radius + j * dx;
            neurons[placedParticles].y = radius + i * dy;
            neurons[placedParticles].radius = radius;
            neurons[placedParticles].synapses = 0;
            placedParticles++;
         }
       }
   }
}


/**************************************************************************************************************
* @brief Function that prints into a file the positions of generated neurons.
* @param file      The file name.
* @param neurons   Array of structures of neurons.
* @param N_neurons Total amount of neurons.
**************************************************************************************************************/
void writeNeuronData(FILE *file, Neuron *neurons, int N_neurons)
{
   if (!file)
   {
      fprintf(stderr, "Error opening file for writing.\n");
      exit(EXIT_FAILURE);
   }

   // Write header
   fprintf(file, "Xlabel   Ylabel   Radius\n");

   // Write neuron data
   for(int i = 0; i < N_neurons; ++i)
   {
      fprintf(file, "%8.4lf %8.4lf %8.2lf\n", neurons[i].x, neurons[i].y, neurons[i].radius);
   }
}


/**************************************************************************************************************
* @brief Computes the Euclidean distance between two neurons.
*
* @param p1 The first neuron.
* @param p2 The second neuron.
*
* @return The Euclidean distance between two neurons.
**************************************************************************************************************/
double Distance(Neuron p1, Neuron p2)
{
   double distance = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
   if(distance < 1e-4)
   {
     return 0.0;
   }
   return distance;
}


/**************************************************************************************************************
* @brief Compares the distance between three neurons and checks if there is an intermediate one.
*
* @param p1 The current neuron.
* @param p2 A possible intermediate neuron.
* @param p3 The current candidate neuron.
*
* @returns True (1) if p2 is between p1 and p3, else, False (0).
**************************************************************************************************************/
bool isBetween(Neuron p1, Neuron p2, Neuron p3)
{
   double distance = Distance(p1, p3);
   if(distance > Distance(p1, p2) && distance > Distance(p2, p3))
   {
     return true;
   }
   return false;
}


/**************************************************************************************************************
* @brief Checks if a new neuron overlaps with any of the existing neurons.
*
* @param neurons   Array of structures of existing neurons.
* @param N_neurons The number of existing neurons.
* @param newNeuron The new neuron to check for overlap.
*
* @return True (1) if there is no overlap, False (0) if there is overlap.
**************************************************************************************************************/
bool isOverlapping(Neuron *neurons, int N_neurons, Neuron newNeuron)
{
   for(int i = 0; i < N_neurons; ++i)
   {
      if(Distance(neurons[i], newNeuron) < 2 * newNeuron.radius)
      {
         return true;
      }
   }
   return false;
}


/**************************************************************************************************************
* @brief Generates a random position inside the specified boundaries.
*
* @param max The upper limit of the box (either width or height).
* @param min The lower limit for valid placement, typically the radius.
*
* @return A random position within the bounds [min, max - min].
**************************************************************************************************************/
float setInside(float max, float min)
{
   float scale = ((float)rand() / RAND_MAX); // Generate a random float between 0 and 1
   return min + scale * (max - 2 * min);     // Scale and shift the value to the desired range
}


/**************************************************************************************************************
* @brief Calculates the probability of forming a connection between two neurons based on distance.
*
* @param dist          The distance between the two neurons.
* @param intDistFactor A factor based on intermediate distances between neurons.
*
* @return A value between 0 and 1, where larger distances reduces the connection probability.
**************************************************************************************************************/
double connectionProbability(double dist, double intDistFactor)
{
   if(dist >= CUTOFF_RADIUS)
   {
      return 0.0;
   }
   return SCALE_FACTOR * (1 - dist/CUTOFF_RADIUS) * exp(- DECAY_FACTOR * (dist/CUTOFF_RADIUS) - INT_DIST_FACTOR * intDistFactor);
}


/**************************************************************************************************************
* @brief Computes the perpendicular distance between a neuron and a pair of neurons.
*
* @param p1 The first neuron, it connects to the second one
* @param p2 The second neuron, it connects to the first one.
* @param p3 The third neuron.
*
* @return The calculated interdistance between the three neurons.
**************************************************************************************************************/
double interDistance(Neuron p1, Neuron p2, Neuron p3)
{
   double lineLength = Distance(p1 , p2);
   if(lineLength < 1e-9)
   {
     return Distance(p1, p3);
   }
   return fabs((p2.y - p1.y) * p3.x - (p2.x - p1.x) * p3.y + p2.x * p1.y - p2.y * p1.x) / lineLength;
}


/**************************************************************************************************************
* @brief Generates a list of nearby neurons for each neuron, based on a cutoff radius.
*
* @param neurons       Array of neurons.
* @param N_neurons     The total number of neurons.
* @param nearbyNeurons A pointer to an array of lists, containing the indices of nearby neurons for each neuron.
* @param nearbyCounts  A pointer to an array storing the number of nearby neurons for each neuron.
**************************************************************************************************************/
void findNearbyParticles(Neuron *neurons, int N_neurons, int ***nearbyNeurons, int **nearbyCounts)
{
   if (N_neurons <= 0)
   {
      fprintf(stderr, "Error: Number of neurons must be positive\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for nearby neurons list and counts
   *nearbyNeurons = (int **)calloc(N_neurons, sizeof(int *));
   *nearbyCounts = (int *)calloc(N_neurons, sizeof(int));

   if (!*nearbyNeurons || !*nearbyCounts)
   {
      fprintf(stderr, "Failed to allocate memory for nearby neurons or counts\n");
      exit(EXIT_FAILURE);
   }

   for(int i = 0; i < N_neurons; ++i)
   {
      (*nearbyCounts)[i] = 0;
      (*nearbyNeurons)[i] = (int *)calloc(N_neurons, sizeof(int)); // allocate max size initially
   }

   // Populate nearby neurons list for each neuron
   for(int i = 0; i < N_neurons; ++i)
   {
      for(int j = 0; j < N_neurons; ++j)
      {
         if(i != j)
         {
            if(Distance(neurons[i], neurons[j]) <= CUTOFF_RADIUS)
            {
               (*nearbyNeurons)[i][(*nearbyCounts)[i]] = j; // Store index of nearby neuron
               (*nearbyCounts)[i]++;
            }
         }
      }
   }
}


/**************************************************************************************************************
* @brief Frees the memory allocated for the nearby neurons lists.
*
* @param nearbyNeurons The list of nearby neurons for each neuron.
* @param N_neurons     The total number of neurons.
* @param nearbyCounts  The array storing the count of nearby neurons for each neuron.
**************************************************************************************************************/
void freeNearbyParticles(int **nearbyNeurons, int N_neurons, int *nearbyCounts)
{
   // Free the nearbyNeurons array if it is not NULL
   if(nearbyNeurons != NULL)
   {
      for(int i = 0; i < N_neurons; ++i)
      {
         free(nearbyNeurons[i]);
      }
      free(nearbyNeurons);
   }

   // Free the nearbyCounts array if it is not NULL
   if(nearbyCounts != NULL)
   {
      free(nearbyCounts);
   }
}



/**************************************************************************************************************
* @brief Allocates memory for a dynamic 2D array with integer elements and initialize them to 0.
*
* @param rows Rows of the array.
* @param cols Columns of the array.
*
* @return The array with memory allocated and all elements initialized to 0.
**************************************************************************************************************/
int** allocate2DIntArray(int rows, int cols)
{
   // Validate parameters
   if(rows <= 0 || cols <= 0)
   {
      fprintf(stderr, "Invalid dimensions: rows and cols must be positive values\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for the array of pointers
   int** array = calloc(rows, sizeof(int*));
   if(!array)
   {
      fprintf(stderr, "Failed to allocate memory for the rows\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for each row and initialize to 0
   for(int i = 0; i < rows; i++)
   {
      array[i] = calloc(cols, sizeof(int)); // Allocate and initialize to 0
      if(!array[i])
      {
         fprintf(stderr, "Failed to allocate memory for row %d\n", i);

         // Free previously allocated memory before exiting
         for(int j = 0; j < i; j++)
         {
            free(array[j]);
         }
         free(array);

         exit(EXIT_FAILURE);
      }
   }

   return array; // Return the allocated 2D array
}


/**************************************************************************************************************
* @brief Frees allocated memory for a dynamic 2D array with integer elements.
*
* @param array The array to be freed (Freedom!).
* @param cols  Columns of the array.
**************************************************************************************************************/
void free2DIntArray(int** array, int rows)
{
   if(!array)
   {
      return; // If the array is NULL, do nothing
   }

   for(int i = 0; i < rows; i++)
   {
      free(array[i]); // Free each row
   }
   free(array); // Free the array of pointers
}


/**************************************************************************************************************
* @brief Allocates memory for a dynamic 2D array with float elements and initialize them to 0.
*
* @param rows Rows of the array.
* @param cols Columns of the array.
*
* @return The array with memory allocated and all elements initialized to 0.
**************************************************************************************************************/
double** allocate2DArray(int rows, int cols)
{
   // Validate parameters
   if(rows <= 0 || cols <= 0)
   {
      fprintf(stderr, "Invalid dimensions: rows and cols must be positive values\n");
      exit(EXIT_FAILURE);
   }
   
   // Allocate memory for the array of pointers
   double** array = calloc(rows, sizeof(double*));
   if(!array)
   {
      fprintf(stderr, "Failed to allocate memory for the rows\n");
      exit(EXIT_FAILURE);
   }

   // Allocate memory for each row and initialize to 0
   for(int i = 0; i < rows; i++)
   {
      array[i] = calloc(cols, sizeof(double)); // Allocate and initialize to 0
      if (!array[i])
      {
         fprintf(stderr, "Failed to allocate memory for row %d\n", i);

         // Free previously allocated memory before exiting
         for(int j = 0; j < i; j++)
         {
            free(array[j]);
         }
         free(array);
         exit(EXIT_FAILURE);
      }
   }

   return array; // Return the allocated 2D array
}



/**************************************************************************************************************
* @brief Frees allocated memory for a dynamic 2D array with float elements.
*
* @param array The array to be freed (Freedom!).
* @param cols  Columns of the array.
**************************************************************************************************************/
void free2DArray(double** array, int rows)
{
   if(!array)
   {
      return; // If the array is NULL, do nothing
   }

   for(int i = 0; i < rows; i++)
   {
      free(array[i]); // Free each row
   }
   free(array); // Free the array of pointers
}



/**************************************************************************************************************
* @brief Allocates memory for a dynamic 3D array and initialize all elements to 0.
*
* @param x First dimension.
* @param y Second direction.
* @param z Third direction.
*
* @return The array with memory allocated and all elements initialized to 0.
**************************************************************************************************************/
double*** allocate3DArray(int x, int y, int z)
{
   // Allocate memory for the array of pointers to 2D arrays
   double ***array = calloc(x, sizeof(double **));
   if (!array)
   {
     fprintf(stderr, "Failed to allocate memory for the first dimension\n");
     exit(EXIT_FAILURE);
   }

   // Allocate memory for each 2D array and initialize to 0
   for(int i = 0; i < x; i++)
   {
      array[i] = calloc(y, sizeof(double *));
      if(!array[i])
      {
         fprintf(stderr, "Failed to allocate memory for the second dimension\n");
         // Free the already allocated memory
         for(int k = 0; k < i; k++)
         {
            free(array[k]);
         }
         free(array);
         exit(EXIT_FAILURE);
      }

      // Allocate memory for each row of the 2D arrays and initialize to 0
      for(int j = 0; j < y; j++)
      {
         array[i][j] = calloc(z, sizeof(double)); // Use calloc to initialize to 0
         if (!array[i][j])
         {
            fprintf(stderr, "Failed to allocate memory for the third dimension\n");
            // Free already allocated memory
            for(int k = 0; k <= i; k++)
            {
                for(int l = 0; l < y; l++)
                {
                    free(array[k][l]);
                }
                free(array[k]);
            }
            free(array);
            exit(EXIT_FAILURE);
         }
      }
   }

   return array; // Return the allocated 3D array
}



/**************************************************************************************************************
* @brief Frees allocated memory for a dynamic 3D array.
*
* @param array The array to be freed (Freedom!).
* @param x     First dimention.
* @param y     Second dimention.
**************************************************************************************************************/
void free3DArray(double ***array, int x, int y)
{
   for(int i = 0; i < x; i++)
   {
      for(int j = 0; j < y; j++)
      {
         free(array[i][j]); // Free each row of the 2D arrays
      }
      free(array[i]); // Free each 2D array
   }
   free(array); // Free the array of pointers to 2D arrays
}


/**************************************************************************************************************
* @brief Initialices neuron positions in random or lattice distribution
*
* @param neurons   Array of structures of neurons.
* @param N_neurons Total amount of neurons.
* @param L_x       Width of the box.
* @param L_y       Height of the box.
* @param radius    Radius of the neurons.
* @param mode      Random (R/r) or Lattice (L/l).
**************************************************************************************************************/
void initializeNeurons(Neuron *neurons, int N_neurons, double L_x, double L_y, double radius, char mode) 
{
   if(mode == 'R' || mode == 'r')
   {
      placeRandom(N_neurons, L_x, L_y, radius, neurons);
   } 
   else if (mode == 'L' || mode == 'l')
   {
      placeLattice(N_neurons, L_x, L_y, radius, neurons);
   } 
   else
   {
      fprintf(stderr, "Invalid mode\n");
      exit(EXIT_FAILURE);
   }
}