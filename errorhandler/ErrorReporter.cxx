
/*

*/

#include "spacommon/spamesh.h"
#include "MeshingDebugger/MeshingDebugger.h"
#include "errorhandler/ErrorReporter.hxx"
#include "errorhandler/Error.hxx"
#include "register/Register.hxx"
#include "interface_api/IErrorHandler.h"

#include "interface_spa/IMeshEdge.h"
#include "interface_spa/IMeshFace.h"

#include "wksp/CErrorContainer.hxx"
#include "wksp/CFilePrint.h"
#include "wksp/CAppWorkSpace.h"

#include <cstdio>


// Be careful with the distinction between ErrorReporter and IErrorHandler call parameter
// ordering...
// ErrorReporter wants entity id first, entity type next...
// IErrorHandler wants entity type first, entity id next...

static char sMessageString[200];
static char sErrorString[200];
static char sEntityType[200];

void ErrorReporter::report_warning(const char scFuncName[], 
                                   int iErrorNum,
                                   int iEntId,
                                   int iEntTyp )
{
    IErrorHandler* poErrorHandler = Register::getErrorHandler();
    decode_error(iErrorNum, sErrorString);
	sprintf(sMessageString, " %s - REPORT_WARNING : %s entity ID = %d \n",scFuncName,sErrorString,iEntId);
    MeshingDebugger::message(sMessageString);

	// Print relevent Warnings (only) in the summary file 
	if ( (static_cast<int>(spa_eHAS_BAD_ELEMENT_DISTORTION)     == iErrorNum) || 
         (static_cast<int>(spa_eHAS_BAD_MESH_QUALITY)           == iErrorNum) ||
         (static_cast<int>(spa_eHAS_BAD_ELEMENT_SKEW_OR_HEIGHT) == iErrorNum) ||
         (static_cast<int>(spa_eHAS_BAD_ELEMENT_FOR_MAP_MESH)   == iErrorNum) ||
         (static_cast<int>(spa_eHAS_BAD_ELEMENT_ASPECT_RATIO)   == iErrorNum) ||
         (static_cast<int>(spa_eHAS_BAD_ELEMENT_COVERAGE)       == iErrorNum) ||
         (static_cast<int>(spa_eHAS_BAD_ELEMENT_COVERAGE_CRUSHED_CYLINDER)       == iErrorNum) ||
         (static_cast<int>(spa_eNON_MANIFOLD_FACET)             == iErrorNum)    )
	{
		CFilePrint::fileOpen( spa_eMESH_SUMMARY );
		CFilePrint::filePrint( spa_eMESH_SUMMARY, sMessageString);
		CFilePrint::fileClose( spa_eMESH_SUMMARY );
	}
    if ( poErrorHandler )
        poErrorHandler->warning( iErrorNum, iEntTyp, iEntId );
    sprintf(sMessageString,"================================================\n");  
    MeshingDebugger::message(sMessageString);
}

void ErrorReporter::report_fatalError( const char scFuncName[],
                                      int iErrorNum,
                                      int iEntId,
                                      int iEntTyp )
{

	sprintf(sMessageString, "%s\n", scFuncName);  
    MeshingDebugger::message( sMessageString);
    decode_error(iErrorNum, sErrorString);

	sprintf(sEntityType, "GEOMETRY_UNKNOWN");
	// Store any fatal error in its container 
	if ( (static_cast<int>(spa_eGEOM_UNKNOWN) != iEntTyp) && 
		 (iEntId > 0) )
	{
		CErrorContainer *pzErrorContainer = CErrorContainer::getInstance();
		CAppWorkSpace* pzAppWrkSpace = CAppWorkSpace::getInstance();
		IMeshFace* pzMFace=NULL;
		IMeshEdge* pzMEdge=NULL;
		if (pzErrorContainer)
		{
			switch(iEntTyp)
			{
			case static_cast<int>(spa_eVOLUME) :
				sprintf(sEntityType, "VOLUME");
				pzErrorContainer->addBodyError(iErrorNum,iEntId);
				break;
			case static_cast<int>(spa_eSURFACE) :
				sprintf(sEntityType, "SURFACE");

				// Check if this is a repair face or not 
				pzMFace=pzAppWrkSpace->getMeshSurface(iEntId);
				if (pzMFace && (false == pzMFace->isRepairFace()) )
				{
					// Add to the error container existing application entities - no repair faces 
					pzErrorContainer->addFaceError(iErrorNum,iEntId);
				}
				break;
			case static_cast<int>(spa_eEDGE) :
				sprintf(sEntityType, "EDGE");

				// Check if this is a repair edge or not 
				pzMEdge=pzAppWrkSpace->getMeshEdge(iEntId);
				if (pzMEdge && (false == pzMEdge->isRepairEdge()) )
				{
					pzErrorContainer->addEdgeError(iErrorNum,iEntId);
				}
				break;
			case static_cast<int>(spa_eVERTEX) :
				sprintf(sEntityType, "VERTEX");
				pzErrorContainer->addVertexError(iErrorNum,iEntId);
				break;
			default:
				break;
			}
		}
	}

	sprintf(sMessageString, "REPORT_FATAL_ERROR - %s : %s %s TAG/ID = %d \n",scFuncName,sErrorString,sEntityType,iEntId);  
    MeshingDebugger::message( sMessageString);
	// Print any fatal error that were reported in the spa summary file 
    CFilePrint::fileOpen( spa_eMESH_SUMMARY );
	CFilePrint::filePrint( spa_eMESH_SUMMARY, sMessageString);
    CFilePrint::fileClose( spa_eMESH_SUMMARY );

    IErrorHandler* poErrorHandler = Register::getErrorHandler(); 
    if ( poErrorHandler )
        poErrorHandler->fatalError( iErrorNum, iEntTyp, iEntId );
    sprintf(sMessageString,"================================================\n");  
    MeshingDebugger::message( sMessageString);
}

int ErrorReporter::report_yesNo(const char scFuncName[],
                                int iErrorNum,
                                int iEntId,
                                int iEntTyp )
{
    IErrorHandler* poErrorHandler = Register::getErrorHandler();
    sprintf(sMessageString, "%s\n", scFuncName);  
    MeshingDebugger::message( sMessageString);
    decode_error(iErrorNum, sErrorString);
    sprintf(sMessageString, "%s\n", sErrorString);  
    MeshingDebugger::message(sMessageString);
    sprintf(sMessageString,"================================================\n");  
    MeshingDebugger::message( sMessageString);
    if ( poErrorHandler )
    {
        return poErrorHandler->yesNo( iErrorNum, iEntTyp, iEntId );
    }
    else
        return 0;
}

void ErrorReporter::decode_error(int iErrorNum, char ErrorMessage[])
{

    switch(iErrorNum)
    {
		// general Error Type
		{
    case((int)spa_eWORKS_FINE):
        sprintf(ErrorMessage,"No problem. Works fine");
        break;
    case((int)spa_eFAILED_TO_MESH):
        sprintf(ErrorMessage,"Failed to mesh");
        break;
    case((int)spa_eMESHER_NOT_IMPLEMENTED):
        sprintf(ErrorMessage,"Mesher Not Implemented");
        break;
    case((int)spa_eUNSUPPORTED_ELEMTYPE):
        sprintf(ErrorMessage,"Unsupported element type");
        break;
    case((int)spa_eMESH_SIZE_TOO_SMALL):
        sprintf(ErrorMessage,"Mesh size is too small");
        break;
    case((int)spa_eMESH_SIZE_TOO_BIG):
        sprintf(ErrorMessage,"Mesh size is too big");
        break;
    case((int)spa_eTOO_MANY_HARD_POINTS):
        sprintf(ErrorMessage,"Too many hardpoints");
        break;
    case((int)spa_eHARDPOINT_TOO_CLOSE_TO_BOUNDARY):
        sprintf(ErrorMessage,"Hardpoint is too close to boundary");
        break;
    case((int)spa_eFAILED_TO_GET_HARD_POINT):
        sprintf(ErrorMessage,"Failed to get hardpoint");
        break;
    case((int)spa_eFAILED_TO_CREATE_NODE):
        sprintf(ErrorMessage,"Failed to create node");
        break;
    case((int)spa_eFAILED_TO_CREATE_MIDNODES):
        sprintf(ErrorMessage,"Failed to create midnodes");
        break;
    case((int)spa_eFAILED_TO_CREATE_ELEM):
        sprintf(ErrorMessage,"Failed to create element");
        break;
    case((int)spa_eFAILED_TO_CREATE_BCKGRD_MESH):
        sprintf(ErrorMessage,"Failed to create background mesh");
        break;
    case((int)spa_eFAILED_TO_FLATTEN_SURFACE):
        sprintf(ErrorMessage,"Failed to flatten surface");
        break;
    case((int)spa_eNODE_NOT_FOUND):
        sprintf(ErrorMessage,"Node not found");
        break;
    case((int)spa_eELEM_NOT_FOUND):
        sprintf(ErrorMessage,"Element not found");
        break;
    case((int)spa_eVERTEX_NOT_FOUND):
        sprintf(ErrorMessage,"Vertex not found");
        break;
    case((int)spa_eEDGE_NOT_FOUND):
        sprintf(ErrorMessage,"Edge not found");
        break;
    case((int)spa_eFACE_NOT_FOUND):
        sprintf(ErrorMessage,"Face not found");
        break;
    case((int)spa_eVOLUME_NOT_FOUND):
        sprintf(ErrorMessage,"Volume not found");
        break;
    case((int)spa_eLOOP_NOT_FOUND):
        sprintf(ErrorMessage,"Loop not found");
        break;
    case((int)spa_eLOOP_NOT_VALID):
        sprintf(ErrorMessage,"Loop not valid");
        break;		
    case((int)spa_eGRADING_NOT_FOUND):
        sprintf(ErrorMessage,"Grading not found");
        break;		
    case((int)spa_eGRADING_TOO_SMALL):
        sprintf(ErrorMessage,"Grading too small");
        break;
    case((int)spa_eGRADING_TOO_BIG):
        sprintf(ErrorMessage,"Grading too big");
        break;
    case((int)spa_eLINEAR_MESH):
        sprintf(ErrorMessage,"Linear Mesh");
        break;
    case((int)spa_eABORT):
        sprintf(ErrorMessage,"Process Aborted");
        break;
    case((int)spa_eNOT_ENOUGH_MEMORY):
        sprintf(ErrorMessage,"Not enough memory");
        break;
    case((int)spa_ePSEUDO_EDGE_NOT_FOUND):
        sprintf(ErrorMessage,"Failed to get pseudo-edge");
        break;
    case((int)spa_eTOO_MANY_PSEUDO_EDGES):
        sprintf(ErrorMessage,"Too many pseudo-edges");
        break;
    case((int)spa_eNO_ACTION):
        sprintf(ErrorMessage,"No action");
        break;
    case((int)spa_eFAILED_TO_MAP_NODE):
        sprintf(ErrorMessage,"Failed to map nodes");
        break;
	}	
		// geometry error Types
		{
    case((int)spa_eLOOP_NOT_CLOSED):
        sprintf(ErrorMessage,"Loop is not closed");
        break;	
    case((int)spa_eFAILED_TO_CREATE_LOOPS):
        sprintf(ErrorMessage,"Failed to create loops");
        break;
    case((int)spa_eFAILED_TO_GET_CURVATURE):
        sprintf(ErrorMessage,"Failed to get curvature");
        break;	
    case((int)spa_eEDGE_TOO_SMALL):
        sprintf(ErrorMessage,"Edge is too small");
        break;
    case((int)spa_eEDGE_TOO_SHORT):
        sprintf(ErrorMessage,"Edge is too short (Collapse short curve) ");
        break;
    case((int)spa_eFACE_TOO_SMALL):
        sprintf(ErrorMessage,"Face is too small");
        break;	
    case((int)spa_eFACE_HAS_NO_EDGE):
        sprintf(ErrorMessage,"Face has no edge");
        break;
    case((int)spa_eFROZEN_EDGE_WITH_NO_NODES ):
        sprintf(ErrorMessage, "Frozen edge has no nodes" );
        break; 
		}

		// 1D Mesher Error types
		{
    case((int)spa_eFAILED_TO_MESH_EDGE):
        sprintf(ErrorMessage,"Failed to mesh edge");
        break;	
    case((int)spa_eFAILED_TO_PROJECT_NODES_ON_EDGE):
        sprintf(ErrorMessage,"Failed to project nodes on edge");
        break;	
    case((int)spa_eFAILED_TO_GET_LOCAL_LENGTH):
        sprintf(ErrorMessage,"Failed to get local length");
        break;	
    case((int)spa_eFAILED_TO_HONOR_LOCAL_LENGTH):
        sprintf(ErrorMessage,"Failed to honor local length");
        break;	
    case((int)spa_eFAILED_TO_REFINE_SAMPLING):
        sprintf(ErrorMessage,"Failed to refine spapling");
        break;	
    case((int)spa_e1DMESHER_NOT_FOUND):
        sprintf(ErrorMessage,"Failed to find 1D mesher");
        break;	
		}

		// 2D Mesher Error types
		{
    case((int)spa_eFAILED_TO_SURFACE_MESH_IN_2D):
        sprintf(ErrorMessage,"Failed to surface mesh in 2D");
        break;	
    case((int)spa_eNO_MESH_DEFINITION):
        sprintf(ErrorMessage,"Failed to create outer loop");
        break;	
    case((int)spa_eFAILED_TO_CREATE_INNER_LOOP):
        sprintf(ErrorMessage,"Failed to create inner loop");
        break;	
    case((int)spa_eFAILED_TO_CREATE_OUTER_LOOP):
        sprintf(ErrorMessage,"Failed to create outer loop");
        break;	
    case((int)spa_eFAILED_TO_GET_WELD_LINE):
        sprintf(ErrorMessage,"Failed to get weld line");
        break;	
    case((int)spa_eFAILED_TO_JOIN_LOOPS):
        sprintf(ErrorMessage,"Failed to join loops");
        break;	
    case((int)spa_eFAILED_TO_CREATE_SIZEMAP):
        sprintf(ErrorMessage,"Failed to compute sizemap");
        break;	
    case((int)spa_eFAILED_TO_PAVE_LOOP):
        sprintf(ErrorMessage,"Failed to pave loop");
        break;	
    case((int)spa_eLOOPS_INTERSECT):
        sprintf(ErrorMessage,"Loops intersect");
        break;	
    case((int)spa_eSELFINTERSECTING_LOOP):
        sprintf(ErrorMessage,"Self-intersecting loops");
        break;
    case((int)spa_eFAILED_TO_GENERATE_SPLIT_LINES):
        sprintf(ErrorMessage,"Failed to generate split lines");
        break;		
    case((int)spa_eFAILED_TO_STORE_QUAD):
        sprintf(ErrorMessage,"Failed to store quad");
        break;		
    case((int)spa_eFAILED_TO_STORE_TRIA):
        sprintf(ErrorMessage,"Failed to store tria");
        break;		
    case((int)spa_eUNREASONABLE_ELEMENT_LENGTH):
        sprintf(ErrorMessage,"Unreasonable element length");
        break;		
    case((int)spa_eELEMENTS_CONNECTED_TO_NODE_NOT_FOUND):
        sprintf(ErrorMessage,"Elements connected to node not found");
        break;		
    case((int)spa_eFAILED_TO_CREATE_FREE_MAPPED_MESH):
        sprintf(ErrorMessage,"Failed to create a mapped-mesh on this surface");
        break;	
    case((int)spa_eBAD_MESH_CONDITION_NUMBER):
        sprintf(ErrorMessage, "Low mesh condition number");
        break;
    case((int)spa_eBAD_MESH_STRETCH):
        sprintf(ErrorMessage, "Low 3D Minimum Stretch (<0.1)");
        break;
    case((int)spa_eFAILED_2D_TO_3D_NODE_MAP):
        sprintf(ErrorMessage, "Failed to map 2D node to 3D space");
        break;
    case((int)spa_eCOINCIDENT_LOOP_NODES ):
        sprintf(ErrorMessage, "Coincident nodes found in 2D loop" );
        break; 
    case((int)spa_eFAILED_TO_CREATE_MAPPED_MESH ):
        sprintf(ErrorMessage, "failed to create the map mesh" );
        break; 
    case((int)spa_eFAILED_STACK_OVERFLOW):
        sprintf(ErrorMessage,"STACK OVERFLOW IN TQM");
        break;
    case((int)spa_eFAILED_TO_MESH_FOR_NO_PARITY):
        sprintf(ErrorMessage,"Failed to mesh for no parity");
        break;
    case((int) spa_eUNMESH_NEXT_TO_NEVER_TINY_EDGE):
        sprintf(ErrorMessage,"Unmesh this face adjacent to a never tiny edge");
        break;
        }

		// 3D Mesher error types.
		{
    case((int)spa_eFAILED_TO_VOLUME_MESH):
        sprintf(ErrorMessage,"Failed to volume mesh");
        break;	
    case((int)spa_eFAILED_TO_SURFACE_MESH_IN_3D):
        sprintf(ErrorMessage,"Failed to surface mesh in 3D");
        break;	
    case((int)spa_eFAILED_TO_TET_MESH):
        sprintf(ErrorMessage,"Failed to Tet mesh");
        break;	
    case((int)spa_eFAILED_TO_HEX_MESH):
        sprintf(ErrorMessage,"Failed to Hex mesh");
        break;	
    case((int)spa_eFAILED_TO_GET_LOOPS):
        sprintf(ErrorMessage,"Failed to get loops");
        break;	
    case((int)spa_eFAILED_TO_GET_FROZEN_DATA):
        sprintf(ErrorMessage,"Failed to get frozen data");
        break;	
    case((int)spa_eFAILED_TO_GET_CONSTRAINT):
        sprintf(ErrorMessage,"Failed to get constraint");
        break;	
    case((int)spa_eFAILED_TO_STORE_MESH):
        sprintf(ErrorMessage,"Failed to store mesh");
        break;	
    case((int)spa_eFAILED_TO_GENERATE_TRIANGLE_MAP):
        sprintf(ErrorMessage,"Failed to generate trianglemap");
        break;	
    case((int)spa_eTRIANGLEMAP_ERROR):
        sprintf(ErrorMessage,"Trianglemap error");
        break;
    case((int)spa_eFAILED_TO_DISCRETIZE_BOUNDARY):
        sprintf(ErrorMessage,"Failed to discretize boundary");
        break;
    case((int)spa_eLENGTH_SCALE_FACTOR_ERROR):
        sprintf(ErrorMessage,"Length scale factor error");
        break;
    case((int)spa_eFAILED_TO_STORE_NODES):
        sprintf(ErrorMessage,"Failed to store nodes");
        break;
    case((int)spa_eELEMENT_NORMAL_ERROR):
        sprintf(ErrorMessage,"Element normal error");
        break;
    case((int)spa_eFAILED_TO_EXTRACT_SURFACE_TRIANGLES):
        sprintf(ErrorMessage,"Failed to extract surface triangles");
        break;
    case((int)spa_eFREE_EDGE_IN_FROZEN_SURFACE_MESH):
        sprintf(ErrorMessage,"Free edge in frozen surface mesh");
        break;
    case((int)spa_eFAILED_INRIA_CHECKS):
        sprintf(ErrorMessage,"Failed INRIA checks");
        break;
    case((int)spa_eFAILED_FROZEN_DATA):
        sprintf(ErrorMessage,"Failed Frozen data");
        break;
    case((int)spa_eFACE_COLLAPSED_BY_FROZEN_DATA):
        sprintf(ErrorMessage,"Face collapsed (zero area) due to Frozen data");
        break;
		}

		// Mesh Cleaner error types.
		{
    case((int)spa_eBAD_INPUT):
        sprintf(ErrorMessage,"Bad input");
        break;	
    case((int)spa_eFAILED_TO_CLEAN_MESH):
        sprintf(ErrorMessage,"Failed to clean mesh");
        break;	
    case((int)spa_eMESH_NOT_CLEANED):
        sprintf(ErrorMessage,"Mesh is not cleaned");
        break;	
    case((int)spa_eMESH_CHANGED):
        sprintf(ErrorMessage,"Mesh has changed");
        break;	
		}

		// Mesh Smoother error types.
		{
    case((int)spa_eFAILED_TO_SMOOTH):
        sprintf(ErrorMessage,"Failed to smooth");
        break;	
    case((int)spa_eINVALID_INPUT_MESH):
        sprintf(ErrorMessage,"Invalid input mesh");
        break;	
    case((int)spa_eELEMENT_EDGE_NOT_FOUND):
        sprintf(ErrorMessage,"Element edge not found");
        break;	
    case((int)spa_eEDGE_SWAP_FAILED):
        sprintf(ErrorMessage,"Edge swap failed");
        break;
		}

		// Facet quality types
		{
    case((int)spa_eGOOD_FACET):
        sprintf(ErrorMessage, "Good facets ");
        break;
    case((int)spa_eNON_MANIFOLD_FACET):
        sprintf(ErrorMessage, "Non-manifold facets present");
        break;
    case((int)spa_eTINY_EDGED_FACET):
        sprintf(ErrorMessage, "Facets with tiny edges present");
        break;
    case((int)spa_eCOLLAPSED_FACET):
        sprintf(ErrorMessage, "Collapsed facets present");
        break;
    case((int)spa_eSELF_INTERSECTING_FACET_LOOPS):
        sprintf(ErrorMessage, "Self-intersecting facet loops present");
        break;
    case((int)spa_eOVERLAPPING_FACET):
        sprintf(ErrorMessage, "Overlapping facets present");
        break;
    case((int)spa_eFAILED_BDRY_DISCR):
        sprintf(ErrorMessage, "Self-intersecting bdry discretisation");
        break;
    case((int)spa_eCOINCIDENT_FACET_NODES):
        sprintf(ErrorMessage,"Coincident Facet nodes found");
        break;
    case((int)spa_eINVALID_LOOP_TOPOLOGY):
        sprintf(ErrorMessage,"Invalid loop topology found");
        break;
		}	
		// Mesh Generation Messages
		{
    case(static_cast<int>(SAM_MESHING_IN_PROGRESS)):
        sprintf(ErrorMessage,"SAM_MESHING_IN_PROGRESS");
        break;
    case(static_cast<int>(SAM_BEGIN_MESHING_PROCESS)):
        sprintf(ErrorMessage,"SAM_BEGIN_MESHING_PROCESS");
        break;
    case(static_cast<int>(SAM_END_MESHING_PROCESS)):
        sprintf(ErrorMessage,"SAM_END_MESHING_PROCESS");
        break;
    case(static_cast<int>(SAM_ABORTING_MESHING_PROCESS)):
        sprintf(ErrorMessage,"SAM_ABORTING_MESHING_PROCESS");
        break;
    case(static_cast<int>(SAM_SUCCESSFULLY_MESHED)):
        sprintf(ErrorMessage,"SAM_SUCCESSFULLY_MESHED");
        break;
    case(static_cast<int>(SAM_LOADING_AND_FLATTENING_VOLUME_FACES)):
        sprintf(ErrorMessage,"SAM_LOADING_AND_FLATTENING_VOLUME_FACES");
        break;
    case(static_cast<int>(SAM_BEGIN_SAMPLE_MESHING)):
        sprintf(ErrorMessage,"SAM_BEGIN_SAMPLE_MESHING");
        break;
    case(static_cast<int>(SAM_END_SAMPLE_MESHING)):
        sprintf(ErrorMessage,"SAM_END_SAMPLE_MESHING");
        break;
    case(static_cast<int>(SAM_BEGIN_SURFACE_MESHING)):
        sprintf(ErrorMessage,"SAM_BEGIN_SURFACE_MESHING");
        break;
    case(static_cast<int>(SAM_END_SURFACE_MESHING)):
        sprintf(ErrorMessage,"SAM_END_SURFACE_MESHING");
        break;
    case(static_cast<int>(SAM_BEGIN_VOLUME_MESHING)):
        sprintf(ErrorMessage,"SAM_BEGIN_VOLUME_MESHING");
        break;
    case(static_cast<int>(SAM_END_VOLUME_MESHING)):
        sprintf(ErrorMessage,"SAM_END_VOLUME_MESHING");
        break;
    case(static_cast<int>(SAM_BEGIN_GENERATING_MIDNODES)):
        sprintf(ErrorMessage,"SAM_BEGIN_GENERATING_MIDNODES");
        break;
    case(static_cast<int>(SAM_DONE_GENERATING_MIDNODES)):
        sprintf(ErrorMessage,"SAM_DONE_GENERATING_MIDNODES");
        break;
    case(static_cast<int>(SAM_SMOOTHING_MESH)):
        sprintf(ErrorMessage,"SAM_SMOOTHING_MESH");
        break;
    case(static_cast<int>(SAM_END_SMOOTHING_MESH)):
        sprintf(ErrorMessage,"SAM_END_SMOOTHING_MESH");
        break;
    case(static_cast<int>(SAM_CLEANING_MESH)):
        sprintf(ErrorMessage,"SAM_CLEANING_MESH");
        break;
    case(static_cast<int>(SAM_END_CLEANING_MESH)):
        sprintf(ErrorMessage,"SAM_END_CLEANING_MESH");
        break;
    case(static_cast<int>(SAM_BEGIN_INTERVAL_ASSIGNMENT)):
        sprintf(ErrorMessage,"SAM_BEGIN_INTERVAL_ASSIGNMENT");
        break;
    case(static_cast<int>(SAM_END_INTERVAL_ASSIGNMENT)):
        sprintf(ErrorMessage,"SAM_END_INTERVAL_ASSIGNMENT");
        break;
		}


		// Tet Meshing Error
		{
    case(static_cast<int>(spa_eINRIA_FAILED_TO_TET_MESH)):
        sprintf(ErrorMessage,"INRIA FAILED TO TET MESH");
        break;
    case(static_cast<int>(spa_eINRIA_OUT_OF_MEMORY)):
        sprintf(ErrorMessage,"INRIA OUT OF MEMORY");
        break;
    case(static_cast<int>(spa_eINRIA_SKIN_INTERSECTS)):
        sprintf(ErrorMessage,"INRIA SKIN INTERSECTS");
        break;
    case(static_cast<int>(spa_eINRIA_BOUNDARY_REGENERATION_FAILED)):
        sprintf(ErrorMessage,"INRIA BOUNDARY REGENERATION FAILED");
        break;
    case(static_cast<int>(spa_eINRIA_DUPLICATE_FACE)):
        sprintf(ErrorMessage,"INRIA DUPLICATE FACE");
        break;
    case(static_cast<int>(spa_eINRIA_PROBABLE_NON_COMPATIBLE_SURF_MESH)):
        sprintf(ErrorMessage,"INRIA PROBABLE NON COMPATIBLE SURF MESH");
        break;
		}

		// Mesh quality status - These are mostly internal errors
		{
    case(static_cast<int>(spa_eMESH_QUALITY_IS_FINE)):
        sprintf(ErrorMessage,"MESH QUALITY IS FINE");
        break;
    case(static_cast<int>(spa_eWAS_GAP_SMOOTHED)):
        sprintf(ErrorMessage,"WAS GAP SMOOTHED");
        break;
	case(static_cast<int>(spa_eHAS_BAD_MESH_CONDITION_NUMBER)):
        sprintf(ErrorMessage," HAS BAD MESH CONDITION NUMBER");
        break;
	case(static_cast<int>(spa_eHAS_BAD_ELEMENT_DISTORTION)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT DISTORTION");
        break;
	case(static_cast<int>(spa_eHAS_BAD_ELEMENT_STRETCH)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT STRETCH");
        break;
	case(static_cast<int>(spa_eHAS_BAD_ELEMENT_WARP)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT WARP");
        break;
	case(static_cast<int>(spa_eHAS_BAD_MESH_QUALITY)):
        sprintf(ErrorMessage,"HAS BAD MESH QUALITY");
        break;
    case(static_cast<int>(spa_eHAS_BAD_ELEMENT_SKEW_OR_HEIGHT)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT SKEW OR HEIGHT");
        break;
    case(static_cast<int>(spa_eHAS_BAD_ELEMENT_FOR_MAP_MESH)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT FOR MAP MESH");
        break;
    case(static_cast<int>(spa_eHAS_BAD_ELEMENT_ASPECT_RATIO)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT ASPECT RATIO");
        break;
    case(static_cast<int>(spa_eHAS_BAD_ELEMENT_COVERAGE)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT COVERAGE");
        break;
    case(static_cast<int>(spa_eHAS_BAD_ELEMENT_COVERAGE_CRUSHED_CYLINDER)):
        sprintf(ErrorMessage,"HAS BAD ELEMENT COVERAGE CRUSHED CYLINDER");
        break;
		}

		// Pyramid errors
		{
    case(static_cast<int>(spa_ePYRAMID_INSERTION_FAILED)):
        sprintf(ErrorMessage,"PYRAMID INSERTION FAILED");
        break;
    case(static_cast<int>(spa_ePYRAMID_OPEN_FAILED)):
        sprintf(ErrorMessage,"PYRAMID OPEN FAILED");
        break;
    case(static_cast<int>(spa_ePYRAMID_MERGE_FAILED)):
        sprintf(ErrorMessage,"PYRAMID MERGE FAILED");
        break;
    case(static_cast<int>(spa_ePYRAMID_SMOOTHING_FAILED)):
        sprintf(ErrorMessage,"PYRAMID SMOOTHING FAILED");
        break;
    case(static_cast<int>(spa_ePYRAMID_SWAPPING_FAILED)):
        sprintf(ErrorMessage,"PYRAMID SWAPPING FAILED");
        break;
    case(static_cast<int>(spa_ePYRAMID_QUALITY_FAILED)):
        sprintf(ErrorMessage,"PYRAMID QUALITY FAILED");
        break;
		}


		// Sweep errors
		{
	case( (int)spa_eSWEEP_SOURCE_TARGET_ARE_IDENTICAL ):
        sprintf(ErrorMessage, "SWEEP SOURCE TARGET ARE IDENTICAL" );
        break;
	case( (int)spa_eSWEEP_SOURCE_TARGET_SHARE_AN_EDGE ):
        sprintf(ErrorMessage, "SWEEP SOURCE TARGET SHARE AN EDGE" );
        break;
	case( (int)spa_eSWEEP_FAILED_TO_MESH_SOURCE_WITH_SCAR ):
        sprintf(ErrorMessage, "SWEEP FAILED TO MESH SOURCE WITH SCAR" );
        break;
	case( (int)spa_eSWEEP_FAILED_TO_MESH_TARGET_WITH_SCAR ):
        sprintf(ErrorMessage, "SWEEP FAILED TO MESH TARGET WITH SCAR" );
        break;
	case( (int)spa_eSWEEP_FAILED_TO_MAP_MESH_WALL_WITH_SCAR):
        sprintf(ErrorMessage, "SWEEP FAILED TO MAP MESH WALL WITH SCAR" );
        break;
    case(static_cast<int>(spa_eSWEEP_NUM_LAYER_CONFLICT)):
        sprintf(ErrorMessage,"NUMBER OF SWEEP LAYER MISMATCH BETWEEN LOOPS");
        break;
    case(static_cast<int>(spa_eSWEEP_NODE_LINE_FAILED)):
        sprintf(ErrorMessage,"SWEEP NODE LINE FAILED");
        break;
    case(static_cast<int>(spa_eSWEEP_FAILED_TO_SMOOTH_LAYER)):
        sprintf(ErrorMessage,"Failed to smooth layer");
        break;
    case(static_cast<int>(spa_eSWEEP_QUALITY_FAILED)):
        sprintf(ErrorMessage,"SWEEP QUALITY FAILED");
        break;
    case(static_cast<int>(spa_eSWEEP_TARGET_PROJECTION_FAILED)):
        sprintf(ErrorMessage,"SWEEP TARGET PROJECTION FAILED");
        break;
    case( (int)spa_eSWEEP_FAILED_TO_MAP_MESH_WALL ):
        sprintf(ErrorMessage, "Failed to map mesh Wall Face during Sweeping" );
        break; 
    case((int)spa_eSWEEP_IGNORED_WALL_HARDPOINTS ):
        sprintf(ErrorMessage, "Ignored Hard Points on Walls during HexSweep" );
        break; 
    case((int)spa_eSWEEP_FAILED_TO_MAP_FROZEN_TARGET ):
        sprintf(ErrorMessage, "SWEEP FAILED TO MAP FROZEN TARGET" );
        break; 
    case((int)spa_eSWEEP_FAILED_EDGE_SEEDING_CONFLICT ):
        sprintf(ErrorMessage, "SWEEP FAILED EDGE SEEDING CONFLICT" );
        break; 
    case((int)spa_eSWEEP_FAILED_TO_MESH_SOURCE ):
        sprintf(ErrorMessage, "SWEEP FAILED TO MESH SOURCE" );
        break; 
    case((int)spa_eSWEEP_FAILED_TO_MESH_TARGET ):
        sprintf(ErrorMessage, "SWEEP FAILED TO MESH TARGET" );
        break; 
    case(static_cast<int>(spa_eSWEEP_FAILED_NUM_LAYERS) ):
        sprintf(ErrorMessage, "FAILED TO COMPUTE NUMBER OF LAYERS " );
        break; 
    case((int)spa_eSWEEP_POOR_SOURCE_MESH_QUALITY ):
        sprintf(ErrorMessage, "Source Face has Poor Quality during HexSweep" );
        break; 
    case((int)spa_eSWEEP_POOR_TARGET_MESH_QUALITY ):
        sprintf(ErrorMessage, "Target Face has Poor Quality during HexSweep" );
        break; 
    case((int)spa_eSWEEP_WALL_MAP_MESH_OVERCONSTRAINED ):
        sprintf(ErrorMessage, "Map mesh wall faces are overconstrained" );
        break; 
		}

		// LocalReMesh errors
		{  
	case((int)spa_eLRM_FAILED_TO_CREATE_SUB_VIRTUAL_EDGES ):
		sprintf(ErrorMessage, "failed to create sub and virtual edges during LocalReMesh" );
        break; 
	case((int)spa_eLRM_FAILED_TO_CREATE_VIRTUAL_LOOPS ):
		sprintf(ErrorMessage, "failed to create virtual loops during LocalReMesh" );
        break;

	case((int)spa_eLRM_FAILED_TO_CREATE_VIRTUAL_LOOPS_DUE_TO_TINY_EDGES ):
		sprintf(ErrorMessage, "failed to create virtual loops due to tiny edges during LocalReMesh" );
        break;

	case((int)spa_eLRM_FAILED_TO_BUILD_VIRTUAL_TOPOLOGY ):
		sprintf(ErrorMessage, "failed to build virtual topology during LocalReMesh" );
        break;

	case((int)spa_eLRM_FAILED_TO_MESH ):
		sprintf(ErrorMessage, "failed to mesh during LocalReMesh" );
        break;
		}

	// Sweep Between errors
		{  
	case(static_cast<int>(spa_eSWP_BTW_FAILED_TO_CREATE_WATER_TIGHT_VOLUME) ):
		sprintf(ErrorMessage, "failed to create water tight volume" );
        break; 
		}
	case(static_cast<int>(spa_eFAILED_TO_HONOR_MAPPED_HOLES)):
		sprintf(ErrorMessage, "failed to honor all mapped holes on face");
		break;

	case(static_cast<int>(spa_eFAILED_TO_HONOR_SOME_MAPPED_HOLES)):
		sprintf(ErrorMessage, "failed to honor some mapped holes on face");
		break;
	case(static_cast<int>(spa_eFAILED_TO_HONOR_WELD_ROWS)):
		sprintf(ErrorMessage, "failed to honor all weld rows on face");
		break;

	case(static_cast<int>(spa_eFAILED_TO_HONOR_SOME_WELD_ROWS)):
		sprintf(ErrorMessage, "failed to honor some weld rows on face");
		break;

	case(static_cast<int>(spa_ePARTIALLY_FILLED_WELD_ROWS)):
		sprintf(ErrorMessage, "filled some weld rows partially due to space crunch");
		break;


	case(static_cast<int>(spa_eFAILED_TO_HONOR_WELD_ROWS_ON_SCARLOOPS)):
		sprintf(ErrorMessage, "cannot honor discontinuous weld rows on scarloops");
		break;
		// default
    default:
        sprintf(ErrorMessage,"Error Message %d Not Implemented",iErrorNum);
        break;
    }

}
