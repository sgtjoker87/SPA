/*
 * CSAMVolumeSupport.cpp - implementation for Volume Support class
 * 
 *
 */

#include "CSAMVolumeSupport.h"

// Initialize the static fields here
CSAMVolumeSupport*   CSAMVolumeSupport::m_pzInstance = 0;

IVolumeSupport* CSAMVolumeSupport::createInstance( sam_TzVolumeSupport& zVolumeSupport )
{
    CSAMVolumeSupport::deleteInstance();

    m_pzInstance = new CSAMVolumeSupport(zVolumeSupport);

    return m_pzInstance;
}
        
IVolumeSupport* CSAMVolumeSupport::getInstance()
{
    return m_pzInstance;
}

void CSAMVolumeSupport::deleteInstance()
{
    if ( m_pzInstance )
        delete m_pzInstance;

    m_pzInstance = 0;
}

CSAMVolumeSupport::CSAMVolumeSupport(sam_TzVolumeSupport& zVolumeSupport)
                                 :m_zVolumeSupport(zVolumeSupport)
{
	return;
}


/**
	 * Query the number of surfaces for the volume
	 * @param int iVolumeId, Volume Id
	 * @return int Number of faces in the volume */
int CSAMVolumeSupport::queryNumFaces(int iVolumeId)
{
	    if ( m_zVolumeSupport.pxQueryNumFaces )
			return m_zVolumeSupport.pxQueryNumFaces( iVolumeId );

		return -1;
}

	/** Query the faces of the volume.
	 * @param iVolumeId, Volume Id
	 * @param aiFaces[], Faces of the volume.
	 * @return Number of faces returned in aiFaces[]. */
int CSAMVolumeSupport::queryFaces(int iVolumeId, int aiFaces[])
{
	if ( m_zVolumeSupport.pxQueryFaces )
		return m_zVolumeSupport.pxQueryFaces( iVolumeId,aiFaces );

		return -1;

}

int CSAMVolumeSupport::queryAppId(int iVolumeId)
{
	if ( m_zVolumeSupport.pxQueryAppId )
		return (m_zVolumeSupport.pxQueryAppId( iVolumeId));

		return (iVolumeId);

}

