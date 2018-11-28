#ifndef MATRIXRECONSTRUCTION_H_
#define MATRIXRECONSTRUCTION_H_

#include <vector>
#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>
#include <mutex>
#include <memory>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

class MatrixReconstructionRenderer : public ManyLightsRenderer{
public:
    MatrixReconstructionRenderer() = delete;
    MatrixReconstructionRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, std::pair<std::uint32_t, std::uint32_t> bucket_size, std::uint32_t light_samples, 
        float min_dist, float step_size_factor, float tolerance, float tau, std::uint32_t max_iterations);
    MatrixReconstructionRenderer(const MatrixReconstructionRenderer& other) = delete;
    MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other);
    MatrixReconstructionRenderer& operator = (const MatrixReconstructionRenderer& other) = delete;
    MatrixReconstructionRenderer& operator = (MatrixReconstructionRenderer&& other);
    ~MatrixReconstructionRenderer();

    bool render(Scene* scene);

    void setCancel(bool cancel){
        std::lock_guard<std::mutex> lock(cancel_lock_);

        cancel_ = cancel;
    }

private:
    std::unique_ptr<ManyLightsClusterer> clusterer_;
    std::pair<std::uint32_t, std::uint32_t> bucket_size_;
    std::uint32_t light_samples_;
    float min_dist_, step_size_factor_, tolerance_, tau_;
    std::uint32_t max_iterations_;
    std::mutex cancel_lock_;
    
    bool cancel_;
};

MTS_NAMESPACE_END

#endif