#ifndef MATRIXSEPARATIONRENDERER_H_
#define MATRIXSEPARATIONRENDERER_H_

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>
#include <mutex>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

enum VISIBILITY{VISIBLE = 0, NOT_VISIBLE, P_VISIBLE, P_NOT_VISIBLE, UNKNOWN};

struct RowSample{
    RowSample(const Point3f& p, const Normal& n, std::uint32_t x, std::uint32_t y, std::uint32_t cols, bool in_scene, Intersection intersection) : 
        position(p), normal(n), image_x(x), image_y(y), col_samples(cols, Spectrum(0.f)), visibility(cols, UNKNOWN), intersected_scene(in_scene), its(intersection){
    }

    Point3f position;
    Normal normal;
    std::uint32_t image_x, image_y;
    std::vector<Spectrum> col_samples;
    std::vector<VISIBILITY> visibility;
    bool intersected_scene;
    Intersection its;
};

class MatrixSeparationRenderer : public ManyLightsRenderer{
public:
    MatrixSeparationRenderer(std::unique_ptr<ManyLightsClusterer> clusterer, 
        float min_dist, float sample_percentage, float error_threshold, float reincorporation_density_threshold,
        std::uint32_t slice_size, std::uint32_t max_prediction_iterations, std::uint32_t max_separation_iterations,
        std::uint32_t show_slices, std::uint32_t only_directsamples, bool separate, bool show_error, bool show_sparse,
        std::uint32_t predictor_mask, bool show_rank);
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
    bool show_slices_, show_only_directsamples_;
    bool cancel_;
    std::mutex cancel_lock_;
    std::vector<RowSample> samples_;
    bool separate_, show_error_, show_sparse_;
    std::uint32_t predictor_mask_;
    bool show_rank_;
};

MTS_NAMESPACE_END

#endif