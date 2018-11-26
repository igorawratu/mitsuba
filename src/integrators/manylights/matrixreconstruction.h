#ifndef MATRIXRECONSTRUCTION_H_
#define MATRIXRECONSTRUCTION_H_

#include <vector>
#include <utility>
#include <eigen3/Eigen/Dense>
#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

class MatrixReconstructionRenderer{
public:
    MatrixReconstructionRenderer() = delete;
    MatrixReconstructionRenderer(std::pair<std::uint32_t, std::uint32_t> bucket_size, std::uint32_t light_samples, 
        float min_dist, float step_size_factor, float tolerance, float tau, std::uint32_t max_iterations);
    MatrixReconstructionRenderer(const MatrixReconstructionRenderer& other);
    MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other);
    MatrixReconstructionRenderer& operator = (const MatrixReconstructionRenderer& other);
    MatrixReconstructionRenderer& operator = (MatrixReconstructionRenderer&& other);
    ~MatrixReconstructionRenderer();

    bool Render(const std::vector<VPL>& vpls, Scene* scene);

    void setCancel(bool cancel){
        cancel_ = cancel;
    }

private:
    std::pair<std::uint32_t, std::uint32_t> bucket_size_;
    std::uint32_t light_samples_;
    float min_dist_, step_size_factor_, tolerance_, tau_;
    std::uint32_t max_iterations_;
    
    bool cancel_;
};

MTS_NAMESPACE_END

#endif