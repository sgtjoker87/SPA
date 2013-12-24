/*

*/

#ifndef FMUTCLEAN_HXX
#define FMUTCLEAN_HXX



#include "cleaner/fmuclean.hxx"
#include "cleaner/nodarray.hxx"
#include "cleaner/elmtarray.hxx"
#include "cleaner/fmuvector.hxx"
#include <cstdio>

class nodClass;
class elmTClass;
class elmQClass;
class fmuEdge;
class fmuMaster;

class TriaCleanTool :
    public CleanTool
{
public:
    virtual ~TriaCleanTool();
    TriaCleanTool();

    int CleanMesh( fmuMaster*, fmuVector& );
      
private:
    int clean_tria( elmTClass* );
    
    int neighbor_quads(elmTClass*, nodClass**, elmQClass**, fmuEdge**,
                           nodClass** );
 
    int double_trias( elmTClass*, elmTClass*, fmuEdge*, nodClass* );

    int two_better_trias( elmTClass*, elmTClass*, fmuEdge*, nodClass*, 
                              nodClass*, nodClass*, nodClass* );

    int tria_and_quad( elmTClass*, elmQClass*, fmuEdge*, nodClass* );

    int two_of_each( elmTClass*, elmQClass*, fmuEdge*, elmTClass*,
                         fmuEdge*, nodClass*, nodClass*, nodClass*,
                         nodClass*, nodClass* );
  
    int nose_to_nose( elmTClass*, elmQClass*, fmuEdge*, elmTClass*,
                          fmuEdge*, nodClass*, nodClass*, nodClass*,
                          nodClass*, nodClass* );
 
    int tria_quad_tria( elmTClass*, elmQClass*, fmuEdge*, elmTClass*,
                            fmuEdge*, nodClass*, nodClass*, nodClass*,
                            nodClass*, nodClass*);
  
    int tria_quad_quad( elmTClass*, elmQClass*, fmuEdge*, elmQClass*,
                            fmuEdge*, nodClass*, nodClass*, nodClass*,
                            nodClass*, nodClass* );
  
    int tria_hourglass( elmTClass*, elmQClass*, fmuEdge*, elmQClass*,
                            fmuEdge*, nodClass*, nodClass*, nodClass*,
                            nodClass*, nodClass*, nodClass*, nodClass* );

    int reshape_quad( elmTClass*, elmQClass*, fmuEdge*, nodClass* );
   
    void fill_6( nodClassDynArray& );

	int check_fill_6( nodClassDynArray& );

    int fill_5( nodClassDynArray&, elmTClass*, elmQClass*, fmuEdge* );

    void delete_this_tria( elmTClass* );

    void delete_this_quad( elmQClass* );
    
    void delete_this_edge( fmuEdge* );

    void delete_this_node( nodClass* );

    void isolate_trias();

    bool check_interior( nodClassDynArray*, int,
                         nodClassDynArray*, int );
   
    void compress_loop( nodClassDynArray* );

	double test_quad_quality(double dx1, double dy1,
		                     double dx2, double dy2,
							 double dx3, double dy3,
							 double dx4, double dy4);

	double get_tria_area (double dx1, double dy1,
                          double dx2, double dy2,
                          double dx3, double dy3);

private:
    elmTClassDynArray   CTriaList;

private:
    TriaCleanTool( const TriaCleanTool& );
    TriaCleanTool& operator= ( const TriaCleanTool& );
};

#endif
