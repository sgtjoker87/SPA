/*
 * CSPAElementSupport.h - class defintion for Element Support
 *
 */


#ifndef CSPAELEMENTSUPPORT_HPP
#define CSPAELEMENTSUPPORT_HPP

#include "interface_api/IElementSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPAElementSupport :
    public IElementSupport 
{
public:
    SPA_DECLARE_NEW_DELETE( CSPAElementSupport )

public:
    static IElementSupport* createInstance( SPA_TzElementSupport *pzElementSupport );
    static IElementSupport* getInstance();
    static void deleteInstance();

public:
    // interface
    virtual bool isStraightenMidNodes( const int shape, const double xyz[][3] ) const;
    virtual bool isSplitQuad( const double xyz[][3] ) const;
    virtual bool isTetraTooFlat( const double xyz[][3] ) const;
    virtual bool doIPassAllEQC( const int shape, const int count, const double xyz[][3] ) const;

public:
    // deprecated
	virtual bool doesElementFailTetCollapse(int       shape,
                                            int       count,
                                            double    xyz[][3],
											double    *pdTetCollapse) const;

	virtual bool doesElementFailAspectRatio(int       shape,
                                            int       count,
                                            double    xyz[][3],
											double    *pdAspectRatio) const;

protected:

    CSPAElementSupport( SPA_TzElementSupport *pzElementSupport );
    virtual ~CSPAElementSupport(){}

private:
    static CSPAElementSupport* m_pzInstance;

private:
	SPA_TzElementSupport* m_pzElementSupport;

private:
    //disable
    CSPAElementSupport( const CSPAElementSupport& );
    CSPAElementSupport& operator=( const CSPAElementSupport& );
};
#endif

