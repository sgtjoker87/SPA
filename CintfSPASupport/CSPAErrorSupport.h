/*
 * CSPAErrorSupport.h - implementation for Error Support class
 *
 */



#ifndef CSPAERRORSUPPORT_HPP
#define CSPAERRORSUPPORT_HPP

#include "interface_api/IErrorHandler.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAErrorSupport : public IErrorHandler
{
public:
    
    SPA_DECLARE_NEW_DELETE( CSPAErrorSupport )

    static IErrorHandler* createInstance(SPA_TzErrorSupport& zErrorSupport);

    static IErrorHandler* getInstance();

    static void deleteInstance();

	/* Fatal Error Call Back */
    void fatalError( int iErrorNum, int iEntType, int iEntId );

    /* Warning Call Back */
    void warning( int iErrorNum, int iEntType, int iEntId );

    /* Prompt (YesNo) Call Back */
    int yesNo ( int iErrorNum, int iEntType, int iEntId );

  
protected:

    CSPAErrorSupport(SPA_TzErrorSupport& zErrorSupport);

    virtual ~CSPAErrorSupport(){};

private:

    static CSPAErrorSupport*    m_pzInstance;

    SPA_TzErrorSupport& m_zErrorSupport;

    // CSPAErrorSupport( const CSPAErrorSupport &X){}

    CSPAErrorSupport& operator= ( 
        const CSPAErrorSupport &X){ return *this; }

};
#endif //CSPAERRORSUPPORT_HPP








