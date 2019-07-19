#ifndef MANYLIGHTSBASE_H_
#define MANYLIGHTSBASE_H_

#include <vector>
#include <mitsuba/render/vpl.h>

MTS_NAMESPACE_BEGIN

class ManyLightsClusterer{
public:
    virtual std::vector<VPL> getClusteringForPoint(const Intersection& its) = 0;
};

class ManyLightsRenderer{
public:
    virtual bool render(Scene* scene, std::uint32_t spp, const RenderJob *job) = 0;
    virtual void setCancel(bool cancel) = 0;
};

MTS_NAMESPACE_END

#endif