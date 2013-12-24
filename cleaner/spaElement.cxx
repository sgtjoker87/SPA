#include "cleaner/spaElement.hxx"

spaElement::spaElement(spaElementType type_, const int *nodes_)
{
	type = type_;
	elementDescriptor = SPA_NORMAL_ELEMENT;
	if (type == SPA_ELEMENT_TRI3)
		numNodes = 3;
	else if (type == SPA_ELEMENT_QUAD4)
		numNodes = 4;
	for (int i=0; i<numNodes; i++)
		nodes[i] = nodes_[i];
	flags = 0x0;
	markValid();
}

void spaElement::getNodes(int &numNodes_, int nodes_[4]) const
{
	numNodes_ = numNodes;
	for (int i=0; i<numNodes; i++)
		nodes_[i] = nodes[i];
}

void spaElement::replaceNode(int from, int to)
{
	for(int i=0; i<numNodes; i++)
	{
		if(nodes[i]==from)
		{
			nodes[i] = to;
		}
	}
}