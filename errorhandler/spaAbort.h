#ifndef SPAABORT_H
#define SPAABORT_H

#include "errorhandler/Error.hxx"

class MeshDriverSPA;

class spaAbort
{
    friend class MeshDriverSPA;

public:

    static void checkUserAbort() throw(Error);

    static void resetUserAbort();

protected:

    static void setUserAbort(bool bAbortFlag);

    static bool m_bAbortFlag;
};

#endif
