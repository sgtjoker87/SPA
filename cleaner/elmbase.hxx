/*
//
//*/

#ifndef ELMBASE_HXX
#define ELMBASE_HXX

#include <SamTypes.h>

#include <cstdio>   // for FILE

class nodClass;
class fmuEdge;

class elmBase
{
public:
    virtual ~elmBase();
    elmBase();

public:
    // elmBase interface
    virtual SAM_ElementShape type() const = 0;

    virtual void nodes( nodClass **apCNodes ) const = 0;
    virtual nodClass* prev_node( nodClass* ) const = 0;
    virtual nodClass* next_node( nodClass* ) const = 0;

    virtual void edges( fmuEdge** ) const = 0;

    virtual void draw_me( int ) const = 0;

public:
    // base class services
    void frozen( bool lNewValue ) { mIsFrozen = lNewValue; }
    bool frozen() const { return mIsFrozen; }

    void id( int iNewId ) { mId = iNewId; }
    int  id() const { return mId; }

    void dump( FILE* ps ) const;

protected:
    int     mId;
    mutable int mindexCurrentNode;
    bool    mIsFrozen;

private:
    // disable
    elmBase( const elmBase& );
    elmBase& operator=( const elmBase& );
};

#endif
