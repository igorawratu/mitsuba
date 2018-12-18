#ifndef MATRIXSEPARATIONRENDERER_H_
#define MATRIXSEPARATIONRENDERER_H_

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>
#include <mutex>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

class MatrixSeparationRenderer : public ManyLightsRenderer{
public:
    MatrixSeparationRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, 
        float min_dist, float sample_percentage, float error_threshold, float reincorporation_density_threshold,
        std::uint32_t slice_size, std::uint32_t max_prediction_iterations, std::uint32_t max_separation_iterations,
        std::uint32_t show_slices);
    MatrixSeparationRenderer(const MatrixSeparationRenderer& other) = delete;
    MatrixSeparationRenderer(MatrixSeparationRenderer&& other);
    MatrixSeparationRenderer& operator = (const MatrixSeparationRenderer& other) = delete;
    MatrixSeparationRenderer& operator = (MatrixSeparationRenderer&& other);
    ~MatrixSeparationRenderer();

    bool render(Scene* scene);

    void setCancel(bool cancel){
        std::lock_guard<std::mutex> lock(cancel_lock_);

        cancel_ = cancel;
    }

private:
    std::unique_ptr<ManyLightsClusterer> clusterer_;
    float min_dist_, sample_percentage_, error_threshold_, reincorporation_density_threshold_;
    std::uint32_t slice_size_, max_prediction_iterations_, max_separation_iterations_;
    bool show_slices_;
    bool cancel_;
    std::mutex cancel_lock_;
};

MTS_NAMESPACE_END

#endif