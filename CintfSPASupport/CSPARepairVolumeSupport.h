/*
 * CSPARepairVolumeSupport.h - class defintion for Repair Volume Support
 *
 *
 */

#ifndef CSPAREPAIRVOLUMESUPPORT_HPP
#define CSPAREPAIRVOLUMESUPPORT_HPP

#include "interface_api/IRepairVolumeSupport.h"
#include "interface_api/SPAInterface.h"
#include "SPAstlmemory/SPANew.hxx"

class CSPARepairVolumeSupport : public IRepairVolumeSupport 
{
public:

    SPA_DECLARE_NEW_DELETE( CSPARepairVolumeSupport )

    static IRepairVolumeSupport* createInstance( SPA_TzRepairVolumeSupport&  zRepairVolumeSupport );
    static IRepairVolumeSupport* getInstance();

    static void deleteInstance();

    virtual int  queryElementSize  ( int iRepVolId, double* pdElemSize );
    virtual int  queryElementType  ( int iRepVolId );
    virtual int  queryElementOrder ( int iRepVolId );
	virtual int  queryNumLayers (int iRepVolId); 
	

protected:

             CSPARepairVolumeSupport( SPA_TzRepairVolumeSupport&  zRepairVolumeSupport );
    virtual ~CSPARepairVolumeSupport(){};

private:

    SPA_TzRepairVolumeSupport& m_zRepairVolumeSupport;

private:

    static CSPARepairVolumeSupport* m_pzInstance;

private:

//           CSPAMeshVolumeSupport( const CSPAMeshVolumeSupport& X ) {}

    CSPARepairVolumeSupport& operator=( const CSPARepairVolumeSupport& X ) { return *this; }

};
#endif //SPAREPAIRVOLUMESUPPORT_HPP

