//====================================================================================================

#include "cleaner/spaMeshImprove.hxx"
#include "utility_math/typeVector.h"
#include "utility_math/mutmathlib.h"
#include "FileImportExport/CSamUnvWrite.h"

char  s_psFileName[256] = { 0 };

#define MESH_IMPROVE_DEBUG 0

double computeEdgeLength(double n0[3], double n1[3])
{
	typeVector3 zV( n0, n1);
	//Vector v(n0, n1);
	//return (typeVector3v.length();
	double dL = zV.Magnitude();
	
	return (dL);
}

spaMeshImprove::spaMeshImprove(spaMesh *m, double minLength)
{
	mesh = m;
	minEdgeLength = minLength;
	faceID = 0;
}

spaMeshImprove::~spaMeshImprove()
{
}

int spaMeshImprove::cleanup()
{
	// Print input mesh info
#if MESH_IMPROVE_DEBUG    
	
    // Write out a universal file of the data
    sprintf(s_psFileName,"spaMeshImproveInput%d.unv",faceID);
    CSamUnvWrite inUnv(s_psFileName, true );
    inUnv.write((*mesh));

	std::cout<<"Input mesh information: "<<std::endl;
	mesh->print();
#endif
	// Mark the boundary nodes
	markBoundaryNodes();

	// Create edges
	createEdges();

	// Collapse edges
	collapseEdges();

#if MESH_IMPROVE_DEBUG  
    // Write out a universal file of the data
    sprintf(s_psFileName,"spaMeshImproveOutput%d.unv",faceID);
    CSamUnvWrite outUnv(s_psFileName, true );
    outUnv.write((*mesh));
#endif

	return 0;
}

spaMesh *spaMeshImprove::getMesh()
{
	return mesh;
}

void spaMeshImprove::markBoundaryNodes()
{
	// Loop over elemenets
	for (int i=0; i<mesh->getNumElement(); i++)
	{
		spaElement &element = mesh->getElement(i);
		// Get nodes
		int numNodes = 0;
		int nodes[4] = {0};
		element.getNodes(numNodes, nodes);
		// Loop over edges
		for (int j=0; j<numNodes; j++)
		{
			// Get an edge (a,b)
			int a = nodes[j];
			int b;
			if (j == (numNodes-1))
				b = nodes[0];
			else
				b = nodes[j+1];

			// Get elements shared by an edge
			std::vector<int> elements;
			mesh->getEdgeNeighbors(a, b, elements);
			if (elements.size() == 1)
			{
				spaNode &nodea = mesh->getNode(a);
				nodea.markBoundary();
				spaNode & nodeb = mesh->getNode(b);
				nodeb.markBoundary();
			}
		}
	}
}

void spaMeshImprove::createEdge(int a, int b, double l)
{
	spaEdge edge(a,b,l);
	edges.push(edge);
}

void spaMeshImprove::createEdges()
{
	// Create small edges that cause element failure due to min edge length
	for (int i=0; i<mesh->getNumElement(); i++)
	{
		spaElement &element = mesh->getElement(i);
		// Get nodes
		int numNode = 0;
		int nodes[4] = {0};
		element.getNodes(numNode, nodes);
		// Loop over edges
		for (int j=0; j<numNode; j++)
		{
			// Get an edge (a,b)
			int a = nodes[j];
			int b;
			if (j == (numNode-1))
				b = nodes[0];
			else
				b = nodes[j+1];

			// Get elements shared by an edge
			std::vector<int> elements;
			mesh->getEdgeNeighbors(a, b, elements);

			// Consider only interior edges
			if (elements.size() == 2)
			{
				// Check if both elements adjacent to an edge are NORMAL to avoid collapsing weld or maped loop elements
				spaElement &element0 = mesh->getElement(elements[0]);
				spaElement &element1 = mesh->getElement(elements[1]);
				if (element0.getElementDiscriptor() == SAM_NORMAL_ELEMENT && element1.getElementDiscriptor() == SAM_NORMAL_ELEMENT)
				{
					// Create only one edge
					if (a < b)
					{
						spaNode &nodea = mesh->getNode(a);
						spaNode &nodeb = mesh->getNode(b);
						// If one boundary node and one interior node 
						if ((nodea.isBoundary() && !nodeb.isBoundary()) 
							|| (!nodea.isBoundary() && nodeb.isBoundary()))
						{
							// If edge length is less than min length, create an edge
							double coord[2][3];
							nodea.getCoord(coord[0]);
							nodeb.getCoord(coord[1]);
							double edgeLength = computeEdgeLength(coord[0], coord[1]);
							if ( edgeLength < 1.15*minEdgeLength)
								createEdge(a, b, edgeLength);
						}
						// If both are interior node - an interior edge
						if (!nodea.isBoundary() && !nodeb.isBoundary())
						{
							// If edge length is less than min length, create an edge
							double coord[2][3];
							nodea.getCoord(coord[0]);
							nodeb.getCoord(coord[1]);
							double edgeLength = computeEdgeLength(coord[0], coord[1]);
							if ( edgeLength < 0.1*minEdgeLength)
								createEdge(a, b, edgeLength);
						}
					}
				}
			}
		}
	}

	// Print edges created
#if MESH_IMPROVE_DEBUG
	std::cout<<edges.size()<<" edges created..."<<std::endl;
#endif
}

bool spaMeshImprove::collapseEdge(int a, int b)
{
	// Find out the elemenets that are going to be deleted and changed
	std::vector<int> deadElements;
	std::vector<int> deltaElements;
	mesh->markNeighbors(b, 0);
	mesh->markNeighbors(a, 1);
	mesh->getMarkedNeighbors(b, deadElements);
	mesh->getUnmarkedNeighbors(b, deltaElements);

	// Make sure edge collapse doesn't create inverted or very small area elements
	const double TOLERANCE = 1e-03;
    for(int i=0; i<deltaElements.size(); i++)
    {
		int ielement = deltaElements[i];
		spaElement &element = mesh->getElement(ielement);

		// Get nodes 
		int numNodes = 0;
		int inodes[4];
		element.getNodes(numNodes, inodes);

		// Virtually replace node b with a
		for (int j=0; j<numNodes; j++)
		{
			if (inodes[j] == b)
				inodes[j] = a;
		}
		
		// Check if element is valid 
		if (element.getType() == SAM_ELEMENT_TRI3)
		{
			double coord[3][3] = {0};
			for (int j=0; j<numNodes; j++)
			{
				spaNode &node = mesh->getNode(inodes[j]);
				node.getCoord(coord[j]);
			}

			// Compute tri area 
			double area = MUTAREA2D(coord[0], coord[1], coord[2]);
			if (area < TOLERANCE)
				return false;
		}
		else if (element.getType() == SAM_ELEMENT_QUAD4)
		{
			double coord[4][3] = {0};
			for (int j=0; j<numNodes; j++)
			{
				spaNode &node = mesh->getNode(inodes[j]);
				node.getCoord(coord[j]);
			}

			// Compute quad area 
			double area = MUTAREA2D(coord[0], coord[1], coord[2]) + MUTAREA2D(coord[2], coord[3], coord[0]);
			if (area < TOLERANCE)
				return false;
		}
    }

    // Delete elements shared by an edge (a,b)
	int numQuad = 0;
    for(int i=0; i<deadElements.size(); i++)
	{
		int ielement = deadElements[i];
		mesh->removeElement(ielement);
		spaElement &element = mesh->getElement(ielement);
		// For quad dominant mesh 
		if (element.getType() == SAM_ELEMENT_QUAD4)
		{
			splitQuad(ielement);
			numQuad++;
		}
	}

	// Re-compute to find new dead elements (tris)
	if (numQuad > 0)
	{
		deadElements.clear();
		deltaElements.clear();
		mesh->markNeighbors(b, 0);
		mesh->markNeighbors(a, 1);
		mesh->getMarkedNeighbors(b, deadElements);
		mesh->getUnmarkedNeighbors(b, deltaElements);
	}

	// Now delete the new dead elements (tris)
	if (numQuad > 0)
	{
		for(int i=0; i<deadElements.size(); i++)
		{
			int element = deadElements[i];
			mesh->removeElement(element);
		}
	}

    // Update the elemenets connected to node b
    for(int i=0; i<deltaElements.size(); i++)
    {
		int ielement = deltaElements[i];
		spaElement &element = mesh->getElement(ielement);

		// Replace node b with node a
		element.replaceNode(b, a);
		
		// Add element to node a's neighbors
		spaNode &node = mesh->getNode(a);
		std::vector<int> &elements = node.getElements();
		elements.push_back(ielement);
    }

    // Delete node b
	mesh->removeNode(b);

	// Collapse edge
#if MESH_IMPROVE_DEBUG
	std::cout<<"Edge "<<a<<" "<<b<<" collapsed..."<<std::endl;
#endif

	return true;
}

void spaMeshImprove::splitQuad(int ielement)
{
	// Get quad nodes
	spaElement &element = mesh->getElement(ielement);
	int numNodes = 0;
	int nodes[4] = {0};
	element.getNodes(numNodes, nodes);

	// Create two tri elements
	int tri[3];
	tri[0]=nodes[0]; tri[1]=nodes[1], tri[2]=nodes[2];
	mesh->addElement(SAM_ELEMENT_TRI3, tri);
	tri[0]=nodes[0]; tri[1]=nodes[2], tri[2]=nodes[3];
	mesh->addElement(SAM_ELEMENT_TRI3, tri);
}

bool spaMeshImprove::collapseEdges()
{
	while (!edges.empty())
	{
		const spaEdge &edge = edges.top();

		spaNode &nodea = mesh->getNode(edge.a);
		spaNode &nodeb = mesh->getNode(edge.b);
		if (nodea.isValid() && nodeb.isValid())
		{
			// Collapse the edge and keep the boundary node
			int a, b;
			if (nodea.isBoundary())
			{
				a = edge.a;
				b = edge.b;
			}
			else
			{
				a = edge.b;
				b = edge.a;
			}
			// Collapse an edge
			collapseEdge(a, b);
		}

		// Remove the edge from queue and delete
		edges.pop();
	}

	return true;
}
