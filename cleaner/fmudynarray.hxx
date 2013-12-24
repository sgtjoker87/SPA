/*

//
// Name
//    fmuDynArray - dynamically sized array of pointers of a specific type.
//
// Description
//    Class fmuDynArray represents a dynamic array of pointers of a specific
//    type. It holds a group of ordered elements. There are two methods of
//    use. One is as an array and elements are accessed through a standard
//    array index. The other is as a loop - beginning follows end - and
//    elements are accessed through relationships to a current position.
//    Internally, storage is implemented as an array, assuring array
//    usage performance. The size of the array is increased as
//    necessary when elements are added.
//
//    Note that this class holds pointers to objects of the specified type.
//
//    To create a specific type of a dynamic array class, two steps
//    are necessary.
//    The first step is to declare the specific dynamic array class
//    via a "declare" statement, as in
//
//       declare(fmuDynArray, <type>);
//
//    The second step is to define the implementation routines for the
//    specific dynamic array class via an "implement" statement, as in
//
//       implement(fmuDynArray, <type>);
//
//    The declare statement is given in a header file that is included
//    by clients of the dynamic array class - the same as with any
//    other class defined within a header file.
//
//    The implement statment is given in a code/implementation file, and
//    compiled into an object file (.o). Note that this file needs to
//    include the header file fmudynarrayi.hxx for the implement
//    statement to work.
//
//    'declare' is a macro that expands into the declaration of a class -
//    specifically a dynamic array class that holds elements of type <type>*.
//    This declared class goes by the name <type>DynArray. For example,
//    to declare a dynamic array of pointers to Vectors, the declare
//    statement would be:
//
//       declare(fmuDynArray, Vector);
//
//    and the declared dynamic array class goes by the name
//
//       VectorDynArray
//
//    Therefore, the name 'VectorDynArray' can be used to create an
//    instance of a Vector dynamic array. That is, 'VectorDynArray'
//    can be used just like any other class name.
//
//
//    The members of fmuDynArray are:
//
//    fmuDynArray();
//      Construct an empty dynamic array.
//
//    fmuDynArray(int n);
//      Construct a dynamic array with initial size n.
//
//    fmuDynArray(int n, <type>* val);
//      Construct a dynamic array with initial size n and initialize
//      all elements to val.
//
//    fmuDynArray(const fmuDynArray& other);
//      Copy constructor - copies the pointers from old to new
//
//    void operator=(const fmuDynArray& other);
//      Assignment operator - copies the pointers from old to new
//
//    <type>*& operator[](int ix);
//      An indexing operator that returns a reference to the i'th element.
//      The array will be lengthed if the index is beyond the end of the
//      array. This reference does not affect the current position.
//
//    void append(<type>* item);
//      Appends the given item to the end of the array. If necessary,
//      the size of the array is increased.
//
//    void insert(<type>* item, int ix);
//      Inserts the given item into the array at the given position. All
//      subsequent elements of the array are moved down. If necessary,
//      the size of the array is increased.
//
//    void remove(int ix);
//      Removes the specified element from the array by moving all
//      subsequent elements up.
//
//    void remove(int ix, int iy);
//      Removes the specified block from the array by moving all
//      subsequent elements up. ix is the index of the first position
//      to be removed. iy is the first position afterward that is kept.
//
//    <type>* get();
//      Returns the element at the current loop position.
//
//    void set(int n);
//      Sets the current loop posiiton to n, where n is an array index.
//
//    int step(int n);
//      Moves the current position by n, which can be positive or negative.
//      If n moves the current position beyond either end of the array, the
//      current position will wrap around to the other end of the array.
//      The value of n defaults to 1. The new current loop position is
//     returned.
//
//    <type>* prev(int n);
//      Returns the element at the nth position before the current position.
//      If this would be beyond the end of the array, the position wraps
//      around the array. The value of n defaults to 1.
//
//    <type>* next(int n);
//      Returns the element at the nth position after the current position.
//      If this would be beyond the end of the array, the position wraps
//      around the array. The value of n defaults to 1.
//
//    int find_near(<type>* item)
//      Find the item in the array first by checking array positions before
//      and after the current postion. Set the current position, and return
//      a value that can be used as an array index.
//
//    int find(<type>* item)
//      Find the item in the array, set the current position, and return
//      a value that can be used as an array index.
//
//    void replace(<type>* item)
//      Put the item into the array at the current position.
//
//    int length() const;
//      Returns the number of elements in the array.
//      Note that this gives the number of elements actually appended or
//      inserted into the array; it does not include any of the "unused"
//      elements.
//
//    void resize(int n);
//      Change the size of the array to 'n' number of elements.
//      This can cause the array to either grow or shrink.
//      after resize( 0 ), empty() will return true.
//
//    void clear();
//      clear all entries out of the array.
//      after clear(), empty() will return true.
//
//    bool empty() const;
//      Return whether the dynamic array is empty (true) or contains
//      entries (false).
//
//    void sort();
//      Sort the array so that the pointers are in ascending order. This
//      does not sort by any field of <type>.
//

*/

#ifndef FMUDYNARRAY_HXX
#define FMUDYNARRAY_HXX

#include "samstlmemory/samSTLAllocator.hxx"
#include <vector>
#include <algorithm>
#include <iterator>

template< typename Type >
class fmuDynArray
{
private:
    typedef std::vector< Type, SAM_STL_ALLOC( Type ) > VEC_Type;

public:
    ~fmuDynArray()
    {
        clear();
    }
    fmuDynArray() :
        mvecType(),
        mCursor( 0 )
    {
    }
    // nn: number of Type to reserve (do not set size!)
    fmuDynArray( size_t nn ) :
        mvecType(),
        mCursor( 0 )
    {
        mvecType.reserve( nn );
    }
    // nn: number of Type to create, initialized to 'val'
    fmuDynArray( size_t nn, const Type& val ) :
        mvecType( nn, val ),
        mCursor( 0 )
    {
    }
    fmuDynArray( const fmuDynArray& COther ) :
        mvecType( COther.mvecType ),
        mCursor( COther.mCursor )
    {
    }
    fmuDynArray& operator=( const fmuDynArray& COther )
    {
        if( &COther != this )
        {
            mvecType = COther.mvecType;
            mCursor = COther.mCursor;
        }
        return *this;
    }
    void absorb( fmuDynArray& COther )
    {
        if( &COther != this )
        {
            clear();
            mvecType.swap( COther.mvecType );
            mCursor = 0;
        }
    }
    Type& operator[]( size_t ix )
    {
        if( ( mvecType.empty() ) ||
            ( ix >= mvecType.size() ) )
            mvecType.resize( ix +1, NULL );
        return mvecType[ ix ];
    }
    Type  operator[]( size_t ix ) const
    {
        if( 0 <= ix  &&  ix < mvecType.size() )
            return mvecType[ ix ];
        return Type();
    }
    void append( Type item )
    {
        mvecType.push_back( item );
    }
    void insert( Type item, int ix )
    {
        if( ix < 0 )
            return;
        if( ix >= (int)mvecType.size() )
        {
            mvecType.resize( ix +1 );
            mvecType[ix] = item;
        }
        else
        {
            //VEC_Type::iterator itr( mvecType.begin() +ix );
            //mvecType.insert( itr, item );
            mvecType.insert( ( mvecType.begin() +ix ), item );
        }
    }
    void remove( size_t ix )
    {
        //if( ix < 0 )
        //    return;
        if( ix > mvecType.size() )
            return;
        //VEC_Type::iterator itr( mvecType.begin() +ix );
        //mvecType.erase( itr );
        mvecType.erase( ( mvecType.begin() +ix ) );
#if 0
        if( ix == mCursor )
            mCursor = 0;
        else if( mvecType.empty() )
            mCursor = 0;
        else if( mCursor >= (int)mvecType.size() )
            mCursor = (int)(mvecType.size() -1);
#else
        if( mvecType.empty() )
            mCursor = 0;
        else if( mCursor >= mvecType.size() )
            mCursor = mvecType.size() -1;
#endif
    }
    void remove( size_t ix, size_t iy )
    {
        //if( ix < 0 )
        //    return;
        //if( iy < 0 )
        //    return;
        const size_t count = mvecType.size();
        //if( ix >= (int)count )
        //    return;
        //if( iy >  (int)count )
        //    return;
        if( ix >= count )
            return;
        if( iy >  count )
            return;
        //VEC_Type::iterator itrB( mvecType.begin() +ix );
        //VEC_Type::iterator itrE( mvecType.begin() +iy );
        //mvecType.erase( itrB, itrE );
        mvecType.erase( ( mvecType.begin() +ix ), ( mvecType.begin() +iy ) );
#if 0
        if( ix <= mCursor && mCursor < iy )
            mCursor = 0;
        else if( mvecType.empty() )
            mCursor = 0;
        else if( mCursor >= (int)mvecType.size() )
            mCursor = (int)(mvecType.size() -1);
#else
        if( mvecType.empty() )
            mCursor = 0;
        else if( mCursor >= mvecType.size() )
            mCursor = mvecType.size() -1;
#endif
    }
    Type get()
    {
        return this->operator[]( mCursor );
    }
    void set( size_t index )
    {
        if( 0 <= index && index < mvecType.size() )
            mCursor = index;
    }
    int step( int nn = 1 )
    {
        mCursor = circularCursor( mCursor +nn );
        return (int)mCursor;
    }
    Type prev( int nn = 1 ) const
    {
        const size_t temp = circularCursor( mCursor -nn );
        return mvecType[temp];
    }
    Type next( int nn = 1 ) const
    {
        const size_t temp = circularCursor( mCursor +nn );
        return mvecType[temp];
    }
    int find_near( Type pItem )
    {
        size_t temp = circularCursor( mCursor +1 );
        if( mvecType[temp] == pItem )
        {
            mCursor = temp;
            return (int)mCursor;
        }
        temp = circularCursor( mCursor -1 );
        if( mvecType[temp] == pItem )
        {
            mCursor = temp;
            return (int)mCursor;
        }
        return find( pItem );
    }
    int find( Type pItem )
    {
        if( mvecType.empty() )
        {
            mCursor = 0;
            return -1;
        }
        const size_t count = mvecType.size();
        if( mCursor < count )
            if( mvecType[mCursor] == pItem )
                return (int)mCursor;
        const size_t cursorSave = mCursor;
        for( mCursor = 0; mCursor < count; ++mCursor )
        {
            if( mvecType[mCursor] == pItem )
                return (int)mCursor;
        }
        mCursor = cursorSave;
        return -1;
    }
    void replace( Type pItem )
    {
        mvecType[mCursor] = pItem;
    }
    void resize( size_t nn )
    {
        mvecType.resize( nn );
        if( mCursor >= mvecType.size() )
            mCursor = 0;
    }
    void clear()
    {
        mvecType.clear();
        mCursor = 0;
    }
    bool empty() const
    {
        return mvecType.empty();
    }
    // not needed
    //void sort()
    //{
    //    std::sort( mvecType.begin(), mvecType.end() );
    //}
    // too many places use 'int' to make this 'size_t' right now.
    int length() const
    {
        return (int)mvecType.size();
    }
private:
#if 0
    int circularCursor( int index ) const
    {
        const int count = (int)mvecType.size();
        if( index < 0 )
            index += count;
        else if( index >= count )
            index -= count;
        return index;
    }
#else
    size_t circularCursor( size_t index ) const
    {
        register const size_t count = mvecType.size();
        // if index is immediately valid, return it
        if( index < count )
            return index;
        // check if index is within 'count' of valid (i.e., 'next' case (1 to count) greater than size)
        register size_t temp = index - count;
        if( temp < count )
            return temp;
        // check if index is within '-count' of valid (i.e., 'prev' case...)
        temp = index + count;
        if( temp < count )
            return temp;
        // invalid index... should never get here.
        // return the invalid index.
        return index;
    }
#endif
private:
    VEC_Type mvecType;
    size_t mCursor;
};

#endif
