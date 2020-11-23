#pragma once

#include "TASSERT.h"
#include "SimpleTypes.h"

class Index3
{
private:
    Index_t m_v1;
    Index_t m_v2;
    Index_t m_v3;
public:
    Index3()
    {
        // leave uninitialized
    }
    Index3(
        Index_t v1,
        Index_t v2,
        Index_t v3)
        : m_v1(v1)
        , m_v2(v2)
        , m_v3(v3)
    {
    }
    Index_t const& v1() const { return m_v1; }
    Index_t const& v2() const { return m_v2; }
    Index_t const& v3() const { return m_v3; }
    void setV1(Index_t n1) { m_v1 = n1; }
    void setV2(Index_t n2) { m_v2 = n2; }
    void setV3(Index_t n3) { m_v3 = n3; }

    // for an alphabetical sort
    bool operator <(Index3 const& other) const
    {
        if (m_v1 < other.m_v1)
            return true;
        else if (m_v1 > other.m_v1)
            return false;
        else
        {
            if (m_v2 < other.m_v2)
                return true;
            else if (m_v2 > other.m_v2)
                return false;
            else
                return m_v3 < other.m_v3;
        }
    }
};

typedef Index3 TriNodes;

