#ifndef SPA_NODE_HXX
#define SPA_NODE_HXX

#include <vector>

// Flags
#define NODE_VALID_FLAG 0x01
#define NODE_BOUNDARY_FLAG 0x02

class spaNode
{
public:
	spaNode() {}
	spaNode(double xyz[2]);
	virtual ~spaNode() {}

    void getCoord(double xyz[2]) const;
    unsigned int isValid() const { return getFlag(NODE_VALID_FLAG); }
    void markValid() { setFlag(NODE_VALID_FLAG); }
    void markInvalid() { unsetFlag(NODE_VALID_FLAG); }
	unsigned int isBoundary() const { return getFlag(NODE_BOUNDARY_FLAG); }
	void markBoundary() { setFlag(NODE_BOUNDARY_FLAG); }
	void markNonboundary() { unsetFlag(NODE_BOUNDARY_FLAG); }
	std::vector<int>& getElements() { return elements; }
	const std::vector<int>& getElements() const { return elements; }

protected:
    unsigned int getFlag(unsigned int flag) const {return flags&flag;}
    void setFlag(unsigned int flag) { flags|=flag; }
    void unsetFlag(unsigned int flag) { flags&= ~flag; }
    unsigned char mark() const { return marks; }
    void mark(unsigned char m) { marks = m; }

private:
	double coord[2];
	unsigned char marks, flags;
	std::vector<int> elements;
};

#endif
