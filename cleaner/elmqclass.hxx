/*

*/



#ifndef ELMQCLASS_HXX
#define ELMQCLASS_HXX

#include "cleaner/elmbase.hxx"

class nodClass;
class fmuEdge;
class fmuMaster;

class elmQClass :
    public elmBase
{
public:
    virtual ~elmQClass();
    elmQClass( nodClass*, nodClass*, nodClass*, nodClass*,
               fmuMaster *master=NULL );
      
public:
    // elmBase interface
    virtual SAM_ElementShape type() const;

    virtual void nodes( nodClass **apCNodeLlist ) const;
    virtual void edges( fmuEdge** ) const;

    virtual nodClass* prev_node( nodClass* ) const;
    virtual nodClass* next_node( nodClass* ) const;

    virtual void draw_me( int ) const;

public:
    // specific to elmQClass
    nodClass* opposite_node( nodClass* ) const;
    fmuEdge*  opposite_edge( fmuEdge* ) const;
      
private:
    bool find_index( nodClass* ) const;

private:
    nodClass *apCNodes[4];
    fmuEdge  *apCEdges[4];

private:
    // disable
    elmQClass();
    elmQClass( const elmQClass& );
    elmQClass& operator=( const elmQClass& );
};

#endif
