// mesh cleaning class description file description block -------
//
// Mesh Cleaning
/*

*/



#ifndef FMUQCLEAN_HXX
#define FMUQCLEAN_HXX


#include "cleaner/fmuclean.hxx"
#include "cleaner/nodarray.hxx"
#include "cleaner/elmqarray.hxx"
#include "cleaner/fmuedgar.hxx"
#include "cleaner/fmu.x"

#include <cstdio>

class nodClass;
class elmQClass;
class fmuEdge;

//Structure used for size functions for speed
struct fmu_TzSize
{
    double *adX;
    double *adY;
    double *adZ;
    double *adLen;
};

class QuadCleanTool :
    public CleanTool
{
public:
    
    QuadCleanTool();

    virtual ~QuadCleanTool();
    
    int CleanMesh( fmuMaster*, fmuVector& ); 

    void AnalyzeMesh( fmuMaster* );

    void EraseElements ();

private:

    int size_clean();

    int size_small_clean( nodClass*, nodClass* );

    int size_small_action( nodClass*, nodClassDynArray&,
                           fmuEdgeDynArray&,
                           elmQClassDynArray&, int );

    int size_big_clean( fmuEdge*, nodClass*, nodClass* );
      
    int valence_clean();
    
    int valence_3_clean( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );

    int valence_4_clean( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
    
    int valence_5_clean( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
  
    int valence_6_clean( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_3434( nodClass*, nodClassDynArray&,
                    fmuEdgeDynArray&, elmQClassDynArray& );
  
    int across_remove_face( nodClass*, nodClassDynArray&,
                            fmuEdgeDynArray&, elmQClassDynArray& );
      
    int double_three( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray&, int );
  
    int triple_three( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );
      
    int squeeze_double_row( nodClass*, nodClassDynArray&,
                            fmuEdgeDynArray&, elmQClassDynArray& );
  
    int clean_4333( nodClass*, nodClassDynArray&,
                    fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_435435( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );
  
    int clean_453453( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );

    int double_triangle( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray&,
                         int*, int* );

    int clean_434435445( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_453443544( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
     
    int clean_434434( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray&, int* );
      
    int clean_443443( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray&, int* );
  
    int fix_double_trans( nodClass*, nodClassDynArray&,
                          fmuEdgeDynArray&, elmQClassDynArray&,
                          int );
      
    int clean_435453( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );
  
    int clean_443544354( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray&,
                         int );
      
    int clean_435443544( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray&,
                         int );
      
    int clean_434543544( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
  
    int clean_434453454( nodClass*, nodClassDynArray&,
                         fmuEdgeDynArray&, elmQClassDynArray& );
  
    int reverse_open_face( nodClass*, nodClassDynArray&,
                           fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_5300003( nodClass*, nodClassDynArray&,
                       fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_55443( nodClass*, nodClassDynArray&,
                     fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_503445( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );
      
    int zip_unzip( nodClass*, nodClassDynArray&, fmuEdgeDynArray&,
                   elmQClassDynArray&, int*, int* );
      
    int clean_554443( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );
      
    int clean_534445( nodClass*, nodClassDynArray&,
                      fmuEdgeDynArray&, elmQClassDynArray& );
      
    int rotate_edge( nodClass*, nodClassDynArray&,
                     fmuEdgeDynArray&, elmQClassDynArray&, int );
      
    int across_open_face( nodClass*, nodClassDynArray&,
                          fmuEdgeDynArray&, elmQClassDynArray& );
     
    int crack_open_middle( nodClass*, nodClassDynArray&,
                           fmuEdgeDynArray&, elmQClassDynArray& );
      
    int boundary_clean();
    
    int bndy_flat_clean( nodClass*, elmQClass*, nodClassDynArray* );

    int bndy_flat_4( nodClassDynArray* );

    int bndy_flat_5( nodClass*, nodClassDynArray&,
                     fmuEdgeDynArray&, elmQClassDynArray& );
      
    int bndy_valence_clean( nodClass*, int, nodClassDynArray&,
                            fmuEdgeDynArray&, elmQClassDynArray& );
  
    int bndy_fix_skew( nodClass*, nodClassDynArray&,
                       elmQClassDynArray&, int );
  
    int bndy_rotate_edge( nodClass*, nodClassDynArray&,
                          fmuEdgeDynArray&, elmQClassDynArray&,
                          int );

    int bndy_543( nodClass*, nodClassDynArray&,
                  fmuEdgeDynArray&, elmQClassDynArray& );

    int bndy_543m( nodClass*, nodClassDynArray&,
                   fmuEdgeDynArray&, elmQClassDynArray& );
   
    int bndy_5043334(nodClass*, nodClassDynArray&,
                     fmuEdgeDynArray&, elmQClassDynArray& );
     
    int bndy_5043(nodClass*, nodClassDynArray&,
                  fmuEdgeDynArray&, elmQClassDynArray&, int );

    int bndy_7(nodClass*, nodClassDynArray&,
               fmuEdgeDynArray&, elmQClassDynArray& );

    int shape_clean();

    int fix_inverted_quad( elmQClass*, nodClass** );

    int collapse_inverted( elmQClass*, nodClass* );

    int split_edge( fmuEdge* );

    int fix_bad_angle( elmQClass*, nodClass* );

    bool edges_intersect( fmuEdge*, fmuEdge* );

    int check_to_squeeze( elmQClass*, nodClass*, nodClass*,
                          nodClass*, nodClass* );

    int valence_2_check( nodClassDynArray& );

    int valence_2_clean( nodClass* );

    int valence_2_hardpt( nodClass* );

    int shape_rotate_next( elmQClass*, nodClass*, nodClass*,
                           nodClass*, nodClass*, bool );

    int shape_544( elmQClass*, nodClass*, nodClass*,
                   nodClass*, nodClass* );

    int shape_rotate_prev( elmQClass*, nodClass*, nodClass*,
                           nodClass*, nodClass*, bool );

    int shape_445( elmQClass*, nodClass*, nodClass*,
                   nodClass*, nodClass* );
 
    int shape_454( elmQClass*, nodClass*, nodClass* );

    int shape_444( elmQClass*, nodClass*, nodClass*,
                   nodClass*, nodClass* );
      
    int shape_444_prev( nodClass*, nodClass*, nodClass*, nodClass*,
                         nodClass*, nodClass*, nodClass*, nodClass*,
                        nodClass*, fmuEdge*, fmuEdge*, fmuEdge*,
                        fmuEdge*, elmQClass*, elmQClass*,
                        elmQClass*, elmQClass* );
   
    int shape_444_next( nodClass*, nodClass*, nodClass*, nodClass*,
                        nodClass*, nodClass*, nodClass*, nodClass*,
                        nodClass*, fmuEdge*, fmuEdge*, fmuEdge*,
                        fmuEdge*, elmQClass*, elmQClass*,
                        elmQClass*, elmQClass* );

    int shape_444_lowpeak( nodClass*, nodClass*, nodClass*, nodClass*,
                           nodClass*, nodClass*, nodClass*, nodClass*,
                           nodClass*, fmuEdge*, fmuEdge*, fmuEdge*,
                           fmuEdge*, elmQClass*, elmQClass*, elmQClass*,
                           elmQClass* );

    int shape_444_left3( nodClass*, nodClass*, nodClass*, nodClass*,
                         nodClass*, nodClass*, fmuEdge*, elmQClass*,
                         elmQClass* );

    int shape_444_right3( nodClass*, nodClass*, nodClass*, nodClass*,
                          nodClass*, nodClass*, fmuEdge*, elmQClass*,
                          elmQClass* );
 
    int shape_434( elmQClass*, nodClass*, nodClass*,
                   nodClass*, nodClass* );

    int shape_smooth_prev( elmQClass*, nodClass*, nodClass* );

    int shape_smooth_next( elmQClass*, nodClass*, nodClass* );

    int shape_last_chance( elmQClass*, nodClass*, nodClass*,
                           nodClass*, nodClass* );

    int fill_choice( nodClassDynArray&, fmuEdgeDynArray&,
                     elmQClassDynArray&, nodClassDynArray& );

    int fill_two( nodClassDynArray&, fmuEdgeDynArray&,
                  elmQClassDynArray&, nodClassDynArray& );

    int fill_three( nodClassDynArray&, fmuEdgeDynArray&,
                    elmQClassDynArray&, nodClassDynArray& );

    int fill_choice_side( nodClassDynArray&, fmuEdgeDynArray&,
                          elmQClassDynArray&, nodClassDynArray& );

    int fill_choice_sidem( nodClassDynArray&, fmuEdgeDynArray&,
                           elmQClassDynArray&, nodClassDynArray& );

    int fill_choice_squeeze( nodClassDynArray&, fmuEdgeDynArray&,
                             elmQClassDynArray&, nodClassDynArray& );
      
    void two_middle_nodes( nodClass*,nodClass*, nodClass**, nodClass** );
      
    void delete_this_quad( elmQClass* );
      
    int delete_this_edge( fmuEdge* );
      
    int delete_this_node( nodClass* );
      
    void reclean( nodClass* );
      
    int neighbors( nodClass*, nodClassDynArray&, fmuEdgeDynArray&, elmQClassDynArray& );
      
    int boundary_neighbors( nodClassDynArray*, nodClassDynArray&,
                            fmuEdgeDynArray&, elmQClassDynArray& );
      
    double angle_on_boundary(nodClassDynArray*);
      
    int get_valence( nodClass* );
    
    bool find_on_boundary( nodClass*, nodClassDynArray**, int*, nodClass** );
      
    void rotate_list( int*, int, int );
      
    void rotate_doubles( double*, int, int );
      
    void print_node( FILE*, nodClass*, const char* );
      
    double edge_ratio( fmuEdge* );
      
    void draw_mesh();
    
    void load_size();

    double get_size( double, double, double );
      
    void unload_size();
      
    void draw_hunting( nodClass*, int );

	double test_quad_quality(double dx1, double dy1,
		                     double dx2, double dy2,
							 double dx3, double dy3,
							 double dx4, double dy4);

    double test_quad_quality (elmQClass *pCQuad);

    void get_two_middle_node_location( nodClass* pCNode1, nodClass* pCNode2,
                                       double adNodeX[], double adNodeY[] );

private:
      
    nodClassDynArray   CCleaningList;
    elmQClassDynArray  CShapeList;

    fmu_TzSize      zSize;
    int             nSize;

    bool            lDrawHunting;
    bool            lDrawSize;

private:
    QuadCleanTool( const QuadCleanTool& );
    QuadCleanTool& operator= ( const QuadCleanTool& );
};

#endif
