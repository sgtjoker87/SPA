
#include "errorhandler/spaAbort.h"
#include "errorhandler/Error.hxx"
#include "interface_api/SamErrorTypes.h"

// Initialize the static data
bool spaAbort::m_bAbortFlag = false;

void spaAbort::checkUserAbort() throw (Error)
{
    if ( m_bAbortFlag )
        throw Error((int)spa_eABORT);

    return;
}

void spaAbort::setUserAbort(bool bAbortFlag)
{
    m_bAbortFlag = bAbortFlag;
    return;
}

void spaAbort::resetUserAbort()
{
    m_bAbortFlag = false;
}

