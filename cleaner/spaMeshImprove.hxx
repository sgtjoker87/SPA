
#ifndef SPA_MESH_IMPROVE_H
#define SPA_MESH_IMPROVE_H

#include <queue>
#include "spaMesh.hxx"

// Dynamic edge 

class spaEdge
{
public:
	spaEdge();
	spaEdge(int n0, int n1, double l) 
	{ a = n0; b = n1; length = l;}
	bool operator<(const spaEdge &edge) const
	{
		return length > edge.length;
	}

	int a, b;
	double length;
};

// Mesh improvement/cleanup 

class spaMeshImprove
{
public:
	spaMeshImprove(spaMesh *m, double minLength);
	virtual ~spaMeshImprove();

	// Main cleanup method
	int cleanup();

	spaMesh *getMesh();
	void setFaceID(int id) { faceID = id; }

protected:
	void markBoundaryNodes();
    void createEdge(int a, int b, double l);
    void createEdges();
	bool collapseEdge(int a, int b);
	void splitQuad(int element);
	bool collapseEdges();

private:
	spaMesh *mesh;
	std::priority_queue<spaEdge> edges;
	double minEdgeLength;
	int faceID;
};

#endif
