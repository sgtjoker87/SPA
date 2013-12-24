#ifndef SPA_ELEMENT_HXX
#define SPA_ELEMENT_HXX

enum spaElementType
{
	SPA_ELEMENT_TRI3,
	SPA_ELEMENT_QUAD4
};

enum spaElementDescriptor
{
	SPA_NORMAL_ELEMENT,
	SPA_WELD_ELEMENT,
	SPA_MAPPED_LOOP_ELEMENT
};

// Flags
#define ELEMENT_VALID_FLAG 0x01

class spaElement
{
public:
	spaElement() {}
	spaElement(spaElementType type, const int *nodes);
	virtual ~spaElement() {}

	spaElementType getType() { return type; }
	void setElementDiscriptor(spaElementDescriptor type) { elementDescriptor = type; }
	spaElementDescriptor getElementDiscriptor() { return elementDescriptor; }
    void getNodes(int &numNodes, int nodes[4]) const;
	void replaceNode(int from, int to);

    unsigned int isValid() const { return getFlag(ELEMENT_VALID_FLAG); }
    void markValid() { setFlag(ELEMENT_VALID_FLAG); }
    void markInvalid() { unsetFlag(ELEMENT_VALID_FLAG); }
    unsigned int getFlag(unsigned int flag) const {return flags&flag;}
    void setFlag(unsigned int flag) { flags|=flag; }
    void unsetFlag(unsigned int flag) { flags&= ~flag; }

	unsigned char mark() const { return marks; }
    void mark(unsigned char m) { marks = m; }

private:
	spaElementType type;
	spaElementDescriptor elementDescriptor;
	int numNodes;
	int nodes[4];
	unsigned char marks, flags;
};

#endif
