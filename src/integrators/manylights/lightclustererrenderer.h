#ifndef LIGHTCLUSTERRENDERER_H_
#define LIGHTCLUSTERRENDERER_H_

#include <vector>
#include <algorithm>
#include <mutex>
#include <mitsuba/render/vpl.h>


#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

class LightClustererRenderer : public ManyLightsRenderer{
public:
    LightClustererRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, float min_dist);

    LightClustererRenderer(const LightClustererRenderer& other) = delete;

    LightClustererRenderer(LightClustererRenderer&& other);

    LightClustererRenderer& operator = (const LightClustererRenderer& other) = delete;

    LightClustererRenderer& operator = (LightClustererRenderer&& other);

    ~LightClustererRenderer();

    bool render(Scene* scene);

    void setCancel(bool cancel);

private:
    float min_dist_;
    std::unique_ptr<ManyLightsClusterer> clusterer_;
    bool cancel_;
    std::mutex cancel_lock_;
};

MTS_NAMESPACE_END

#endif