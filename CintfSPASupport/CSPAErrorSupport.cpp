/*
 * CSPAErrorSupport.cpp - implementation for Error Support class
 *
 */


#include "CSPAErrorSupport.h"

// Initialize the static fields here
CSPAErrorSupport*   CSPAErrorSupport::m_pzInstance = 0;

IErrorHandler* CSPAErrorSupport::createInstance( SPA_TzErrorSupport& zErrorSupport )
{
    CSPAErrorSupport::deleteInstance();

    m_pzInstance = new CSPAErrorSupport(zErrorSupport);

    return m_pzInstance;
}
        
IErrorHandler* CSPAErrorSupport::getInstance()
{
    return m_pzInstance;
}

void CSPAErrorSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSPAErrorSupport::CSPAErrorSupport(SPA_TzErrorSupport& zErrorSupport)
                                 :m_zErrorSupport(zErrorSupport)
{
	return;
}


/* Fatal Error Call Back */
void CSPAErrorSupport::fatalError( int iErrorNum, int iEntType, int iEntId )
{
	if (m_zErrorSupport.pxFatalError)
	{
		m_zErrorSupport.pxFatalError(iErrorNum,iEntType,iEntId );
		return;
	}
	else
	{
		return;
	}
}

/* Warning Call Back */
void CSPAErrorSupport::warning( int iErrorNum, int iEntType, int iEntId )
{
	if (m_zErrorSupport.pxWarning)
	{ 
		m_zErrorSupport.pxWarning(iErrorNum,iEntType,iEntId );
		return;
	}
	else
	{
		return;
	}
}

/* Prompt (YesNo) Call Back */
int CSPAErrorSupport::yesNo( int iErrorNum, int iEntType, int iEntId )
{
	if (m_zErrorSupport.pxYesNo)
	{
        return m_zErrorSupport.pxYesNo(iErrorNum,iEntType,iEntId );
	}
	else
	{
		return-1;
	}
}
