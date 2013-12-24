/*
 * CSPAMemorySupport.cpp - implementation for Memory Support
 *
 *
 *
 */

#include "CSPAMemorySupport.h"

// Initialize the static fields here
CSPAMemorySupport*   CSPAMemorySupport::m_pzInstance = 0;

IMemoryAllocator* CSPAMemorySupport::createInstance( SPA_TzMemorySupport& zMemorySupport )
{
    CSPAMemorySupport::deleteInstance();

    m_pzInstance = new CSPAMemorySupport(zMemorySupport);

    return m_pzInstance;
}
        
IMemoryAllocator* CSPAMemorySupport::getInstance()
{
    return m_pzInstance;
}

void CSPAMemorySupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAMemorySupport::CSPAMemorySupport(SPA_TzMemorySupport &zMemorySupport )
                                    : m_zMemorySupport( zMemorySupport )
{
}

void* CSPAMemorySupport::SPA_alloc(size_t n )
{
    if ( m_zMemorySupport.pxSPA_alloc )
        return m_zMemorySupport.pxSPA_alloc(n);

    return NULL;
}

void CSPAMemorySupport::SPA_free(void *P)
{
    if ( m_zMemorySupport.pxSPA_free )
        m_zMemorySupport.pxSPA_free(P);
}


