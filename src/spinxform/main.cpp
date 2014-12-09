// =============================================================================
// SpinXForm -- main.cpp
// Keenan Crane
// August 16, 2011
//

#include <iostream>
#include <fstream>
#include "Mesh.h"
#include "Image.h"

using namespace std;

int main( int argc, char **argv )
{
   if( argc != 4 )
   {
      cerr << "usage: " << argv[0] << " mesh.obj image.tga result.obj" << endl;
      return 1;
   }

   // load mesh
   Mesh mesh;
   mesh.read( argv[1] );

   // load image
   Image image;
   image.read( argv[2] );

   // apply transformation
   const double scale = 5.;
   mesh.setCurvatureChange( image, scale );
   mesh.updateDeformation();

   // write result
   mesh.write( argv[3] );

   return 0;
}


int load_point_curvature(const char *filename,
                         vector<double> &point_curvature)
{
  ifstream ifs(filename);
  if(ifs.fail()){
      cerr << "# [error] can not open point_curvature file." << endl;
      return __LINE__;
    }
  size_t num;
  ifs >> num ;
  point_curvature.resize(num);
  for(size_t i = 0; i  < num; ++i){
      ifs >> point_curvature[i];
    }
  return 0;
}

//int main( int argc, char **argv )
//{
//  if( argc != 4 )
//    {
//      cerr << "usage: " << argv[0] << " mesh.obj point_curvature result.obj" << endl;
//      return 1;
//    }

//  // load mesh
//  Mesh mesh;
//  mesh.read( argv[1] );

//  vector<double> point_curvature;
//  load_point_curvature(argv[2], point_curvature);

//  // apply transformation
//  mesh.setCurvatureChange(point_curvature,5);
//  mesh.updateDeformation();

//  // write result
//  mesh.write( argv[3] );

//  return 0;
//}

