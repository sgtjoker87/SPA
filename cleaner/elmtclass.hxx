/*

*/


#ifndef ELMTCLASS_HXX
#define ELMTCLASS_HXX

#include "cleaner/elmbase.hxx"

class nodClass;
class fmuEdge;
class fmuMaster;

class elmTClass :
    public elmBase
{
public:
    virtual ~elmTClass();
    elmTClass( nodClass*, nodClass*, nodClass*, fmuMaster *master=NULL );
    
public:
    // elmBase interface
    virtual SAM_ElementShape type() const;

    virtual void nodes( nodClass **apCNodeLlist ) const;
    virtual nodClass* prev_node( nodClass* ) const;
    virtual nodClass* next_node( nodClass* ) const;

    virtual void edges(fmuEdge**) const;

    virtual void draw_me( int ) const;

public:
    // specific to elmTClass
    void check_freeze();

private:
    bool find_index( nodClass* ) const;

private:
    nodClass *apCNodes[3];
    fmuEdge  *apCEdges[3];

private:
    // disable
    elmTClass();
    elmTClass( const elmTClass& );
    elmTClass& operator=( const elmTClass& );
};

#endif
