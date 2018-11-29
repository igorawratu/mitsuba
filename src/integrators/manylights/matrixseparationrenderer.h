#ifndef MATRIXSEPARATIONRENDERER_H_
#define MATRIXSEPARATIONRENDERER_H_

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>
#include <mutex>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

class MatrixSeparationRenderer : public ManyLightsRenderer{
public:
    MatrixSeparationRenderer(std::uint32_t slice_size, float min_dist);
    MatrixSeparationRenderer(const MatrixSeparationRenderer& other);
    MatrixSeparationRenderer(MatrixSeparationRenderer&& other);
    MatrixSeparationRenderer& operator = (const MatrixSeparationRenderer& other);
    MatrixSeparationRenderer& operator = (MatrixSeparationRenderer&& other);
    ~MatrixSeparationRenderer();

    bool Render(Scene* scene);

    void setCancel(bool cancel){
        std::lock_guard<std::mutex> lock(cancel_lock_);

        cancel_ = cancel;
    }

private:
    std::uint32_t slice_size_;
    float min_dist_;
    bool cancel_;
    std::mutex cancel_lock_;
};

MTS_NAMESPACE_END

#endif