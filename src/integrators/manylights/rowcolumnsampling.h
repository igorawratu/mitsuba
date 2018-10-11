#ifndef ROWCOLUMNSAMPLING_H_
#define ROWCOLUMNSAMPLING_H_

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

#include <vector>
#include <tuple>

class RowColumnSampling{
public:
    RowColumnSampling() = delete;
    RowColumnSampling(const std::vector<VPL>& vpls, std::uint32_t num_clusters, std::uint32_t rows, std::uint32_t cols,
        std::tuple<std::uint32_t, std::uint32_t> resolution);
    RowColumnSampling(const RowColumnSampling& other);
    RowColumnSampling(RowColumnSampling&& other);
    RowColumnSampling& operator = (const RowColumnSampling& other);
    RowColumnSampling& operator = (RowColumnSampling&& other);
    ~RowColumnSampling();

    void updateClustering(const std::vector<VPL>& vpls, std::uint32_t num_clusters, std::uint32_t rows, std::uint32_t cols,
        std::tuple<std::uint32_t, std::uint32_t> resolution);

    std::vector<VPL> getClusteringForPoint();

private:
    std::vector<VPL> clustering_;
};

#endif