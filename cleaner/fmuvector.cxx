/*

*/

#include "cleaner/fmuvector.hxx"
#include "cleaner/nodclass.hxx"
#include <cmath>

#define  TWOPI    6.28318530717958647692


fmuVector::fmuVector()
{
  
   dVecX = 0.0;
   dVecY = 0.0;
   dVecZ = 0.0;
   
}

fmuVector::fmuVector( nodClass* pCNode )
{
   
   dVecX = pCNode->node_x();
   dVecY = pCNode->node_y();
   dVecZ = pCNode->node_z();
  
}


// Description:
//
// Create a mathematical vector from one node to another
//
// Access:
//
// fmuVector( pCNodeA, pCNodeB );
//
// Input Parameters:
//
// pCNodeA       nodClass*       node coordinates of "from" node
// pCNodeB       nodClass*       node coordinates of "to" node
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
fmuVector::fmuVector( nodClass* pCNodeA, nodClass* pCNodeB )
{
  
   dVecX = pCNodeB->node_x() - pCNodeA->node_x();
   dVecY = pCNodeB->node_y() - pCNodeA->node_y();
   dVecZ = pCNodeB->node_z() - pCNodeA->node_z();
  
}


// Description:
//
// Create a mathematical vector from 3 double precision values
//
// Access:
//
// fmuVector( dNewx, dNewy, dNewz );
//
// Input Parameters:
//
// dNewx          double         Original X value for vector
// dNewy          double         Original Y value for vector
// dNewz          double         Original Z value for vector
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
fmuVector::fmuVector( double dNewx, double dNewy,
                      double dNewz )
{
   
    dVecX = dNewx;
    dVecY = dNewy;
    dVecZ = dNewz;
  
}


// Description:
//
// Copy constructor, create one vector from another
//
// Access:
//
// fmuVector( COther );
//
// Input Parameters:
//
// COther         fmuVector&          Vector to copy from
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
fmuVector::fmuVector( const fmuVector& COther )
{

    dVecX = COther.dVecX;
    dVecY = COther.dVecY;
    dVecZ = COther.dVecZ;
  
}


// Description:
//
// destructor
//
// Access:
//
// delete vector
//
// Input Parameters:
//
// none
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
fmuVector::~fmuVector()
{
   
}


// Description:
//
// Unitize the vector, make it of unit length
//
// Access:
//
// unitize()
//
// Input Parameters:
//
// none
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
void fmuVector::unitize()
{
   
   double dSquare, dMag;
   dSquare = dVecX * dVecX +  dVecY * dVecY +  dVecZ * dVecZ;
   dMag = sqrt( dSquare );
   dVecX /= dMag;
   dVecY /= dMag;
   dVecZ /= dMag;
  
}


// Description:
//
// Compute the squared length of a vector
//
// Access:
//
// len_sq
//
// Input Parameters:
//
// none
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// Additional comments:
// Using the squared length can save time if simply determining whether
// one vector is longer than another.
//
// ------------------------------------------
//
double fmuVector::len_sq() const
{
  
   double dSquare;
   dSquare = dVecX * dVecX +  dVecY * dVecY +  dVecZ * dVecZ;
   return dSquare;
  
}


// Description:
//
// Compute the length of a vector
//
// Access:
//
// length()
//
// Input Parameters:
//
// none
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
double fmuVector::length() const
{
  
   double dSquare;
   dSquare = dVecX * dVecX +  dVecY * dVecY +  dVecZ * dVecZ;
   return sqrt( dSquare );
}


// Description:
//
// Assignment operator, assign the value of one vector to another
//
// Access:
//
// this = COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector with the desired values
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
void fmuVector::operator=( const fmuVector& COther )
{
  
   dVecX = COther.dVecX;
   dVecY = COther.dVecY;
   dVecZ = COther.dVecZ;
  
}


// Description:
//
// Add two vectors together. 
//
// Access:
//
// this + COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector to be added to this one
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The vector containing the sum is returned
//
// ------------------------------------------
//
fmuVector fmuVector::operator+( const fmuVector& COther ) const
{
   
   fmuVector CTemp;
   CTemp.dVecX = dVecX + COther.dVecX;
   CTemp.dVecY = dVecY + COther.dVecY;
   CTemp.dVecZ = dVecZ + COther.dVecZ;

   return CTemp;
   
}


// Description:
//
// Subtract a vector from this vector
//
// Access:
//
// this - COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector to subtract from this one
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The vector containing the difference is returned
//
// ------------------------------------------
//
fmuVector fmuVector::operator-( const fmuVector& COther ) const
{
 
   fmuVector CTemp;
   CTemp.dVecX = dVecX - COther.dVecX;
   CTemp.dVecY = dVecY - COther.dVecY;
   CTemp.dVecZ = dVecZ - COther.dVecZ;

   return CTemp;
  
}


// Description:
//
// Negate a vector
//
// Access:
//
// -this
//
// Input Parameters:
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The vector containing the negation is returned
//
// ------------------------------------------
//
fmuVector fmuVector::operator-() const
{
  
   fmuVector CTemp;
   CTemp.dVecX = -dVecX;
   CTemp.dVecY = -dVecY;
   CTemp.dVecZ = -dVecZ;

   return CTemp;
  
}


// Description:
//
// Multiply this vector by a scalar
//
// Access:
//
// this * dValue
//
// Input Parameters:
//
// dValue         double         the scalar value to multiply by
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The vector containing the scaled values is returned
//
// ------------------------------------------
//
fmuVector fmuVector::operator*( const double dValue ) const
{
  
   fmuVector CTemp;
   CTemp.dVecX = dVecX * dValue;
   CTemp.dVecY = dVecY * dValue;
   CTemp.dVecZ = dVecZ * dValue;

   return CTemp;
  
}


// Description:
//
// Divide this vector by a scalar
//
// Access:
//
// this / dValue
//
// Input Parameters:
//
// dValue         double         the scalar value to divide by
//
// Output Parameters:
//
// none
//
// Return Code:
//
// the vector containing the scaled values is returned
//
// ------------------------------------------
//
fmuVector fmuVector::operator/( const double dValue ) const
{
  
   fmuVector CTemp;
   CTemp.dVecX = dVecX / dValue;
   CTemp.dVecY = dVecY / dValue;
   CTemp.dVecZ = dVecZ / dValue;

   return CTemp;
  
}


// Description:
//
// Add a vector to this vector
//
// Access:
//
// this += COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector to be added to this one
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
void fmuVector::operator+=( const fmuVector& COther )
{
  
   dVecX += COther.dVecX;
   dVecY += COther.dVecY;
   dVecZ += COther.dVecZ;
   
}


// Description:
//
// Subtract a vector from this vector
//
// Access:
//
// this -= COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector to subtract from this one
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
void fmuVector::operator-=( const fmuVector& COther )
{
  
   dVecX -= COther.dVecX;
   dVecY -= COther.dVecY;
   dVecZ -= COther.dVecZ;
  
}


// Description:
//
// Multiply this vector by a scalar
//
// Access:
//
// this *= dValue
//
// Input Parameters:
//
// dValue         double         the scalar value to multiply by
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
void fmuVector::operator*=( const double dValue )
{
  
   dVecX *= dValue;
   dVecY *= dValue;
   dVecZ *= dValue;
  
}


// Description:
//
// Divide this vector by a scalar
//
// Access:
//
// this /= dValue
//
// Input Parameters:
//
// dValue         double         the scalar value to divide by
//
// Output Parameters:
//
// none
//
// Return Code:
//
// none
//
// ------------------------------------------
//
void fmuVector::operator/=( const double dValue )
{
  
   dVecX /= dValue;
   dVecY /= dValue;
   dVecZ /= dValue;
  
}


// Description:
//
// Compute the cross product of two vectors
//
// Access:
//
// this * COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector to cross with this one
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The vector containing the cross product is returned
//
// ------------------------------------------
//
fmuVector fmuVector::operator*( const fmuVector& COther ) const
{
   
   fmuVector CTemp;
   CTemp.dVecX = dVecY * COther.dVecZ - dVecZ * COther.dVecY;
   CTemp.dVecY = dVecZ * COther.dVecX - dVecX * COther.dVecZ;
   CTemp.dVecZ = dVecX * COther.dVecY - dVecY * COther.dVecX;

   return CTemp;
  
}


// Description:
//
// Compute the dot product of two vectors
//
// Access:
//
// this % COther
//
// Input Parameters:
//
// COther         fmuVector&         The other vector to dot with this one
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The double precision value of the cosine between the vectors
//
// ------------------------------------------
//
double fmuVector::operator%( const fmuVector& COther ) const
{
  
   double dDot;
   dDot = dVecX * COther.dVecX + dVecY * COther.dVecY + dVecZ * COther.dVecZ;
   return dDot;
}


// Description:
//
// Compute the angle between two vectors as the counteclockwise
// angle from the first vector to the second. This routine has been
// optimized for speed. It assumes the vectors are in 2D (z component
// is zero). (this*) is assumed to be the normal vector of either
// (0, 0, 1) or (0, 0, -1).
//
// Access:
//
// fmuVector( node );
//
// Input Parameters:
//
// CToNext         fmuVector&         First vector or vector to the next
//                                    node in a loop
// CToPrev         fmuVector&         Second vector or vector to the previous
//                                    node in a loop
//
// Output Parameters:
//
// none
//
// Return Code:
//
// The angle in radians is returned
//
// ------------------------------------------
//
double fmuVector::angle( const fmuVector& CToNext, const fmuVector& CToPrev ) const
{
   double dSineABC=0.0;
   fmuVector CPerp( -(CToNext.dVecY), CToNext.dVecX, 0.0 );
   double dXpart = CToPrev.dVecX * CToNext.dVecX +
                   CToPrev.dVecY * CToNext.dVecY;
   double dYpart = CToPrev.dVecX * CPerp.dVecX +
                   CToPrev.dVecY * CPerp.dVecY;
   if( dXpart == 0.0 && dYpart == 0.0)
        return 0.0;
   double dAngle = atan2( dYpart, dXpart );
   //atan2 has shown roundoff problems. Trap for them so the result
   //isn't 2PI instead of zero.
   if( dAngle > -fmu_ANGLE_TOL && dAngle < fmu_ANGLE_TOL ) 
   {
	   // Check for concavity
	   dSineABC = -CToPrev.dVecX*CToNext.dVecY + CToPrev.dVecY*CToNext.dVecX;
            
	   // If concave, set the angle to 360 deg
	   if (dSineABC <= 1.0e-12)
			dAngle += TWOPI;
	   else
		   return 0.0;
   }
   if( dAngle < 0.0 )
        dAngle += TWOPI;
   if( dVecZ == -1.0)
        dAngle = TWOPI - dAngle;

   return dAngle;
}


// Description:
//
// Rotate a vector through an angle with this vector as the plane
// normal. If the normal is in the -Z direction, use the negative
// of the angle.
//
// Access:
//
// rotate( CVector, dAngle )
//
// Input Parameters:
//
// CVector         fmuVector&         the vector to be rotated
// dAngle          double   the angle in radians
//
// Output Parameters:
//
// none
//
// Return Code:
//
// the rotated vector is returned
//
// ------------------------------------------
//
fmuVector fmuVector::rotate( const fmuVector& CVector,
                             const double dAngle )
{
    fmuVector COutVec;

    double dLocAngle;
    if( dVecZ < 0 )
        dLocAngle = -1.0 * dAngle;
    else
        dLocAngle = dAngle;

    const double dCosine = cos( dLocAngle );
    const double dSine = sin( dLocAngle );

    COutVec.dVecX = ( CVector.dVecX * dCosine ) - ( CVector.dVecY * dSine );
    COutVec.dVecY = ( CVector.dVecX * dSine ) + ( CVector.dVecY * dCosine );
    COutVec.dVecZ = 0.0;

    return COutVec;
}
