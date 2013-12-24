#include "Error.hxx"

Error::Error(int iErrorNum )
            :m_iErrorNum(iErrorNum),
             m_iEntTyp(0),
             m_iEntId(0)
{
}

Error::Error(int iErrorNum,
             int iEntTyp,
             int iEntId )
            :m_iErrorNum(iErrorNum),
             m_iEntTyp(iEntTyp),
             m_iEntId(iEntId)
{
}

Error::~Error()
{
}

int Error::getErrorNum()
{
    return m_iErrorNum;
}

int Error::getEntTyp()
{
    return m_iEntTyp;
}

int Error::getEntId()
{
    return m_iEntId;
}


