/*
 * CSPAMemorySupport.h - class defintion for Memory Support
 *
 *
 */

/* Generated by Together */

#ifndef CSPAMEMORYSUPPORT_HPP
#define CSPAMEMORYSUPPORT_HPP

#include "interface_api/IMemoryAllocator.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAMemorySupport : public IMemoryAllocator 
{
public:

    // SPA_DECLARE_NEW_DELETE( CSPAMemorySupport )

    static IMemoryAllocator* createInstance(SPA_TzMemorySupport&  zMemorySupport);

    static IMemoryAllocator* getInstance();

    static void deleteInstance();

    virtual void* SPA_alloc( size_t n );

    virtual void SPA_free( void  *P );

protected:

    CSPAMemorySupport(SPA_TzMemorySupport&  zMemorySupport);

    virtual ~CSPAMemorySupport(){}

private:

    static  CSPAMemorySupport*  m_pzInstance;

	SPA_TzMemorySupport&    m_zMemorySupport;

    // CSPAFaceSupport( const CSPAFaceSupport &X){}

    CSPAMemorySupport& operator= (
        const CSPAMemorySupport &X){ return *this;}

};
#endif //MEMORYSUPPORTSPA_HPP

