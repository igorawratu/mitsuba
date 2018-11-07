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
    MatrixReconstructionRenderer();
    MatrixReconstructionRenderer(const MatrixReconstructionRenderer& other);
    MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other);
    MatrixReconstructionRenderer& operator = (const MatrixReconstructionRenderer& other);
    MatrixReconstructionRenderer& operator = (MatrixReconstructionRenderer&& other);
    ~MatrixReconstructionRenderer();

    bool Render(Scene* scene, const std::vector<VPL>& vpls, 
        const std::pair<std::uint32_t, std::uint32_t>& bucket_size, const std::uint32_t& light_samples, 
        float min_dist, std::uint8_t* output_image);

    void setCancel(bool cancel){
        cancel_ = cancel;
    }

private:
    bool cancel_;
};

MTS_NAMESPACE_END

#endif