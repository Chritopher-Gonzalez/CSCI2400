#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "Filter.h"

using namespace std;

#include "rdtsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

class Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
	int value;
	input >> value;
	filter -> set(i,j,value);
      }
    }
    return filter;
  } else {
    cerr << "Bad input in readFilter:" << filename << endl;
    exit(-1);
  }
}


double
applyFilter(class Filter *filter, cs1300bmp *input, cs1300bmp *output)
{
    long long cycStart, cycStop;
    cycStart = rdtscll();
    
    output -> width = input -> width;
    output -> height = input -> height;
    /*local variables*/
    
    //Bounds on for loops that tranverse the grid
    int inputWidth = (input -> width) - 1;
    int inputHeight = (input -> height) - 1;
    //Moving getDivisor calls outside of loops
    int filterDivisor = filter -> getDivisor();
    //Moving getSize calls outside of loops
    int filterSize = filter -> getSize();
    
    //loop to transverse grids
    //change the order of first two loops
    // make it faster sometime
    //#pragma omp parallel for
    for(int row = 1; row < inputHeight; row = row + 1) {
        for(int col = 1; col < inputWidth; col = col + 1) {
            int outputRedPlane = 0;
            int outputGreenPlane = 0;
            int outputBluePlane = 0;
            
            // loops to apply filter
            //changed order of filter loops
            for (int j = 0; j < filterSize; j++) {
                //temp variable to save calcuation time
                int c = col + j - 1;
                for (int i = 0; i < filterSize; i++) {
                    //temp variable to save calcuation time
                    int filterData = filter -> get(i, j);
                    int r = row + i - 1;
                    
                    outputRedPlane += (input -> color[0][r][c] * filterData);
                    outputGreenPlane += (input -> color[1][r][c] * filterData);
                    outputBluePlane += (input -> color[2][r][c] * filterData);
                }
            }
            
            outputRedPlane /= filterDivisor;
            outputGreenPlane /= filterDivisor;
            outputBluePlane /= filterDivisor;
            
            //red plane
            outputRedPlane = outputRedPlane < 0 ? 0 : (outputRedPlane > 255 ? 255 : outputRedPlane);
            //green plane
            outputGreenPlane = outputGreenPlane < 0 ? 0 : (outputGreenPlane > 255 ? 255 : outputGreenPlane);
            //blue plane
            outputBluePlane = outputBluePlane < 0 ? 0 : (outputBluePlane > 255 ? 255 : outputBluePlane);
            
            output -> color[0][row][col] = outputRedPlane;
            output -> color[1][row][col] = outputGreenPlane;
            output -> color[2][row][col] = outputBluePlane;
        }
    }
   
    cycStop = rdtscll();
    double diff = cycStop - cycStart;
    double diffPerPixel = diff / (output -> width * output -> height);
    fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n", diff, diff / (output -> width * output -> height));
    return diffPerPixel;
}


                /*                     
                    //temp variable to save calcuation time
                    int filterData1 = filter -> get(0, j);
                    int filterData2 = filter -> get(1, j);
                    int filterData3 = filter -> get(2, j);
                    int r1 = row + 0 - 1;
                    int r2 = row + 1 - 1;
                    int r3 = row + 2 - 1;
                    
                    outputRedPlane += (input -> color[0][r1][c] * filterData1);
                    outputRedPlane += (input -> color[0][r2][c] * filterData2);
                    outputRedPlane += (input -> color[0][r3][c] * filterData3);
                    outputGreenPlane += (input -> color[1][r1][c] * filterData1);
                    outputGreenPlane += (input -> color[1][r2][c] * filterData2);
                    outputGreenPlane += (input -> color[1][r3][c] * filterData3);
                    outputBluePlane += (input -> color[2][r1][c] * filterData1);
                    outputBluePlane += (input -> color[2][r2][c] * filterData2);
                    outputBluePlane += (input -> color[2][r3][c] * filterData3);
                    

                 */
