// =============================================================================
// SpinXForm -- Image.h
// Keenan Crane
// August 16, 2011
//
// Image represents a grayscale bitmapped image.  Standard usage might
// look something like
//
//    Image im;
//    im.load( "image.tga" );
//
//    // sample image at point p
//    double p[2] = { .5, 1.23 };
//    double value = im.sample( p[0], p[1] );
//

#ifndef SPINXFORM_IMAGE_H
#define SPINXFORM_IMAGE_H

#include <vector>
#include <string>

using namespace std;

class Image
{
   public:
            double& operator()( int x, int y );
      const double& operator()( int x, int y ) const;
      // accesses double (x,y)

      double sample( double x, double y ) const;
      // samples image at (x,y) using bilinear filtering

      int  width( void ) const;
      int height( void ) const;
      // returns image dimensions

      void read( const char* filename );
      // loads an image file in Truevision TGA format
      // (must be uncompressed RGB image with 24 or 32 bits per double)

      void reload( void );
      // updates image from disk

   protected:
      void clamp( int& x, int& y ) const;
      // clamps coordinates to range [0,w-1] x [0,h-1]

      string filename; // name of source file
      vector<double> pixels; // interleaved RGBA double data in range [0-1]
      int w; // width
      int h; // height
};

#endif
