#ifndef ERROR_H
#define ERROR_H

#include "spastlmemory/spaNew.hxx"

class Error
{
public:

    SPA_DECLARE_NEW_DELETE( Error )

    Error( int iErrorNum );
    Error( int iErrorNum, int iEntTyp, int iEntId );

    virtual ~Error();

    int getErrorNum();
    int getEntTyp();
    int getEntId();

private:

    int m_iErrorNum;
    int m_iEntTyp;
    int m_iEntId;

  private:
  // private to preclude usage
  // Error    ( const Error& X );
    Error& operator=( const Error& X );
};

#endif // ERROR_H
