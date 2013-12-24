/*
*/

#ifndef NODCLASS_HXX
#define NODCLASS_HXX


#include "cleaner/nodarray.hxx"
#include "cleaner/fmuedgar.hxx"
#include <cstdio>
#include <cstddef>

class fmuEdge;
class fmuVector;
class fmuMaster;

class nodClass
{
public:
    virtual ~nodClass();

    nodClass();
    nodClass( double dnewX, double dnewY, double dnewZ,
              fmuMaster *master=NULL,    
              double dnewU=-1.0,
              double dnewV=-1.0  );
    nodClass( fmuVector&, fmuMaster *pCMaster=NULL,
              double dnewU=-1.0,
              double dnewV=-1.0 );
   
    void set( fmuVector& );
    void set(double dX, double dY, double dZ){dNodeX = dX; dNodeY = dY; dNodeZ = dZ;}
    void set_uv(double u, double v) { dNodeU = u; dNodeV = v;}
    void lock_node() { lLockNode = true; }
    void unlock_node() { lLockNode = false; }
    double node_x() const { return dNodeX;}
    double node_y() const { return dNodeY;}
    double node_z() const { return dNodeZ;}      

    double node_u() const { return dNodeU;}
    double node_v() const { return dNodeV;}

    void hardpt( bool lNewValue ) { lHardpoint = lNewValue; }
    bool hardpt() const { return lHardpoint; }

    void smooth_hardpt( bool lNewValue ) {lSmoothHP  = lNewValue; }
    bool smooth_hardpt() const { return lSmoothHP; }

    void boundary( bool lNewValue ) { lBoundary = lNewValue; }
    bool boundary() const { return lBoundary; }

    void type( int iNewValue ) { miType = iNewValue; }
    int  type() const { return miType; }

    void   gradient( double dNewValue ) { dGradient = dNewValue; }
    double gradient() const { return dGradient; }

    int number_edges() const;   
    bool add_edge( fmuEdge* );
    void delete_edge( fmuEdge* );
    void edge_list( fmuEdgeDynArray& ) const;
    fmuEdge* shared_edge( nodClass* ) const;
    fmuEdge* first_edge() const;

    bool smooth_me( fmuVector* );

    bool equal_angles( fmuEdge*, fmuVector*, fmuVector* );

    void draw_me( int ) const;
    void move_to_me() const;
    void draw_to_me( int ) const;

    double size();

    void id( int iNewId ) { mId = iNewId; }
    int  id() const { return mId; }

    bool in_smooth_list() const { return lSmoothList; }
    void in_smooth_list( bool lFlag) { lSmoothList = lFlag; }

    bool in_clean_list() const { return lCleanList; }
    void in_clean_list( bool lFlag) { lCleanList = lFlag; }

    bool promoted_to_boundary() const { return lPromoted; }
    void promoted_to_boundary( bool lFlag) { lPromoted = lFlag; }

    bool is_on_degenerate_loop() const { return lDegenerateLoopPoint; }
    void is_on_degenerate_loop( bool lFlag) { lDegenerateLoopPoint = lFlag; }      

    void set_proj(nodClass *projNode) { pCProjectNode = projNode;}
    nodClass* get_proj() const {return  pCProjectNode;}

    void surround_nodes(nodClassDynArray &nodeList);

    void set_source( nodClass *sourceNode) {pCSourceNode = sourceNode;}
    nodClass* get_source() const { return pCSourceNode;}

    void dump( FILE* ps ) const;

private:
      
    fmuEdgeDynArray edgeList;
    
    nodClass *pCProjectNode;
    nodClass *pCSourceNode;

    double dGradient;
    double dNodeX;
    double dNodeY;
    double dNodeZ;
    double dNodeU;
    double dNodeV;

    int mId;  
    int miType;

    bool lPromoted;
    bool lDegenerateLoopPoint;
    bool lCleanList;
    bool lSmoothList;
    bool lSmoothHP;
    bool lHardpoint;
    bool lBoundary;
    bool lLockNode;

private:
    //disable
    nodClass( const nodClass& );
    nodClass& operator= ( const nodClass&  );
};

#endif
