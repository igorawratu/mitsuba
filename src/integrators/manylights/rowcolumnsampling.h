#ifndef ROWCOLUMNSAMPLING_H_
#define ROWCOLUMNSAMPLING_H_

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

#include <vector>
#include <tuple>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

class RowColumnSampling : public ManyLightsClusterer{
public:
    RowColumnSampling() = delete;
    RowColumnSampling(const std::vector<VPL>& vpls, std::uint32_t rows, std::uint32_t cols,
        std::tuple<std::uint32_t, std::uint32_t> resolution, const Scene* scene, float min_dist, bool vsl);
    RowColumnSampling(const RowColumnSampling& other);
    RowColumnSampling(RowColumnSampling&& other);
    RowColumnSampling& operator = (const RowColumnSampling& other);
    RowColumnSampling& operator = (RowColumnSampling&& other);
    ~RowColumnSampling();

    std::vector<VPL> getClusteringForPoint(const Intersection& its);

private:
    std::vector<VPL> clustering_;
    float min_dist_;
};

MTS_NAMESPACE_END

#endif