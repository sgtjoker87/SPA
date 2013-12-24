#include "cleaner/spaMesh.hxx"
#include <iostream>
#include <cassert>
#include <algorithm>

spaMesh::spaMesh(int numNodes, int numTris, int numQuads)
{
	nodes.reserve(numNodes);
	elements.reserve(numTris+numQuads);
}

int spaMesh::addNode(double coord[2])
{
	spaNode node(coord);
	nodes.push_back(node);
	int id = (int)nodes.size()-1;
	return id;
}

void spaMesh::removeNode(int n)
{
	spaNode &node = getNode(n);
	node.markInvalid();
	std::vector<int> &elements = node.getElements();
	elements.clear();
}

int spaMesh::addElement(spaElementType type, const int *n)
{
	spaElement element(type, n);
	elements.push_back(element);
	int id = (int)elements.size()-1;
	
	int numNodes = 0;
	if (type == SAM_ELEMENT_TRI3)
		numNodes = 3;
	else if (type == SAM_ELEMENT_QUAD4)
		numNodes = 4;

	for (int i=0; i<numNodes; i++)
	{
		spaNode &node = getNode(n[i]);
		std::vector<int> &elems = node.getElements();
		elems.push_back(id);
	}

	return id;
}

void spaMesh::removeElement(int e)
{
	spaElement &elem = getElement(e);
	elem.markInvalid();

	// Get nodes 
	int numNodes = 0;
	int nodes[4] = {0};
	elem.getNodes(numNodes, nodes);
	int found=0;
    std::vector<int>::iterator itr;
	for (int i=0; i<numNodes; i++)
	{
		// Remove element from node
		spaNode &node = getNode(nodes[i]);
		std::vector<int> &elems = node.getElements();
		itr = find(elems.begin(), elems.end(), e);
		if (itr != elems.end())
		{
			found++; 
			elems.erase(itr);
		}
		// Assert
		assert( found > 0 );
		itr = find(elems.begin(), elems.end(), e);
		assert(itr == elems.end());
	}
}

void spaMesh::getEdgeNeighbors(int n1, int n2, std::vector<int> &outElements)
{
    markNeighbors(n1, 1);
    markNeighbors(n2, 0);
    getUnmarkedNeighbors(n1, outElements);
}

void spaMesh::markNeighbors(int n, unsigned short mark)
{
	spaNode &node = getNode(n);
	std::vector<int> &elems = node.getElements();
    for (int i=0; i<elems.size(); i++) 
    {
		spaElement &elem = getElement(elems[i]);
		elem.mark(mark);
    }
}

void spaMesh::getUnmarkedNeighbors(int n, std::vector<int> &outElements)
{
	spaNode &spanode = getNode(n);
	std::vector<int> &elems = spanode.getElements();
    for (int i=0; i<elems.size(); i++) 
    {
        spaElement &spaelem = getElement(elems[i]);
		if (!spaelem.mark())
		{
			outElements.push_back(elems[i]);
		}
    }
}

void spaMesh::getMarkedNeighbors(int n, std::vector<int> &outElements)
{
	spaNode &spanode = getNode(n);
	std::vector<int> &elems = spanode.getElements();
    for (int i=0; i<elems.size(); i++) 
    {
        spaElement &spaelem = getElement(elems[i]);
		if (spaelem.mark())
		{
			outElements.push_back(elems[i]);
		}
    }
}

void spaMesh::print()
{
	std::cout<<"Numbe of nodes = "<<getNumNode()<<std::endl;
	std::cout<<"Number of elements = "<<getNumElement()<<std::endl;
}