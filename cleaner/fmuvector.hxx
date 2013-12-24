/*

*/


#ifndef FMUVECTOR_HXX
#define FMUVECTOR_HXX

const double fmu_ANGLE_TOL = 0.000000001;

class nodClass;

class fmuVector
{
public:

    virtual ~fmuVector();
      
    fmuVector();
    fmuVector( nodClass* );
    fmuVector( nodClass*, nodClass* );
    fmuVector( double, double, double );
    fmuVector( const fmuVector& );

    void x( double dNewx ) { dVecX = dNewx;}
    double x() const { return dVecX; }

    void y( double dNewy ) { dVecY = dNewy;}
    double y() const { return dVecY; }

    void z( double dNewz ) { dVecZ = dNewz;}
    double z() const { return dVecZ; }

    void unitize();                             // make unit length
    double len_sq() const;                      // compute length squared
    double length() const;                      // compute length

    void operator=( const fmuVector& );                 // assignment operator

    fmuVector operator+( const fmuVector& ) const;      // add two vectors
    fmuVector operator-( const fmuVector& ) const;      // subtract vector
    fmuVector operator-() const;                        // negate vector
    fmuVector operator*( const double ) const;          // multiply by scalar
    fmuVector operator/( const double ) const;          // divide by scalar

    void operator+=( const fmuVector& );               // add. vector to this
    void operator-=( const fmuVector& );               // sub. vector from this
    void operator*=( const double );         // mult. this by scalar
    void operator/=( const double );         // div. this by scalar

    fmuVector operator*( const fmuVector& ) const;    // cross product
    double operator%( const fmuVector& ) const;       // dot product

    fmuVector rotate( const fmuVector&, const double );         // rotate through angle

    // assumes (*this) is the normal to the two vectors.
    double angle( const fmuVector&, const fmuVector& ) const;   // angle between

private:
    double dVecX;
    double dVecY;
    double dVecZ;
};

#endif
