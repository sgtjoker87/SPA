#ifndef SAM_MESH_HXX
#define SAM_MESH_HXX

#include "cleaner/spaNode.hxx"
#include "cleaner/spaElement.hxx"

class spaMesh
{
public:
	spaMesh(int numNodes, int numTris, int numQuads);
	virtual ~spaMesh() {}

    int addNode(double coord[2]);
	void removeNode(int n);
	spaNode& getNode(int n) { return nodes[n]; }
	int getNumNode() const { return (int) nodes.size(); }

    int addElement(spaElementType type, const int *nodes);
	void removeElement(int e);
	spaElement& getElement(int e) { return elements[e]; }
	int getNumElement() const { return (int) elements.size(); }

	void getEdgeNeighbors(int n1, int n2, std::vector<int> &elements);
    void markNeighbors(int n, unsigned short mark=0);
	void getUnmarkedNeighbors(int n, std::vector<int> &elements);
	void getMarkedNeighbors(int n, std::vector<int> &elements);

	void print();

private:
	std::vector<spaNode> nodes;
	std::vector<spaElement> elements;
};

#endif

