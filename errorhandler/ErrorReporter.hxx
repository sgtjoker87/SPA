#ifndef ERRORREPORTER_H
#define ERRORREPORTER_H

#include "spastlmemory/spaNew.hxx"



class MeshDriverSPA;

class ErrorReporter
{
    friend class MeshDriverSPA;

public:

    SPA_DECLARE_NEW_DELETE( ErrorReporter )

    static void report_warning( const char scFuncName[], 
                                int iErrorNum,
                                int iEntId, 
                                int iEntTyp
                              );

    static void report_fatalError( const char scFuncName[],
                                   int iErrorNum,
                                   int iEntId, 
                                   int iEntTyp
                                 );

    static int  report_yesNo( const char scFuncName[],
                              int iErroNum,
                              int iEntId,
                              int iEntTyp
                             );

    static void decode_error( int iErrorNum, char ErrorMessage[] );

};

#endif // ERRORREPORTER_H
