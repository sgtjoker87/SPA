#include "cleaner/spaNode.hxx"

spaNode::spaNode(double xyz[2])
{
    for (int i=0; i<2; i++)
        coord[i] = xyz[i];

	flags = 0x0;
	markValid();
}

void spaNode::getCoord(double xyz[2]) const
{
    for (int i=0; i<2; i++)
        xyz[i] = coord[i];
}
