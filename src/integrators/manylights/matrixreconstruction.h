#ifndef MATRIXRECONSTRUCTION_H_
#define MATRIXRECONSTRUCTION_H_

#include <vector>
#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>
#include <mutex>
#include <memory>

#include "manylightsbase.h"
#include "common.h"

MTS_NAMESPACE_BEGIN

struct ReconstructionSample{
    ReconstructionSample() : intersected_scene(false){
    }

    ReconstructionSample(std::uint32_t x, std::uint32_t y, bool in_scene, Intersection intersection, const Ray& r) : 
        image_x(x), image_y(y), intersected_scene(in_scene), its(intersection), ray(r), color(0.f), 
        fully_sampled_color(0.f){
    }

    ReconstructionSample(ReconstructionSample&& other) : image_x(other.image_x), image_y(other.image_y), 
        intersected_scene(other.intersected_scene), its(other.its), ray(other.ray), color(other.color),
        fully_sampled_color(other.fully_sampled_color), unoccluded_samples(other.unoccluded_samples){
    }

    ReconstructionSample& operator = (ReconstructionSample&& other){
        if(this != &other){
            image_x = other.image_x;
            image_y = other.image_y;
            intersected_scene = other.intersected_scene;
            its = other.its;
            ray = other.ray;
            color = other.color;
            fully_sampled_color = other.fully_sampled_color;
            unoccluded_samples = other.unoccluded_samples;
        }

        return *this;
    }

    std::uint32_t image_x, image_y;
    bool intersected_scene;
    Intersection its;
    Ray ray;
    Spectrum color;
    Spectrum fully_sampled_color;
    std::uint32_t visibility_errors;
    std::vector<Spectrum> unoccluded_samples;
};

enum ClusteringStrategy{MDLC, LS};

class MatrixReconstructionRenderer : public ManyLightsRenderer{
public:
    MatrixReconstructionRenderer() = delete;
    MatrixReconstructionRenderer(const std::vector<VPL>& vpls, float sample_percentage_, 
        float min_dist, float step_size_factor, float tolerance, float tau, std::uint32_t max_iterations,
        std::uint32_t slice_size, bool visibility_only, bool adaptive_col, bool adaptive_importance_sampling, 
        bool adaptive_force_resample, bool adaptive_recover_transpose, bool truncated, bool show_slices, bool vsl,
        bool gather_stat_images, bool show_svd, ClusteringStrategy clustering_strategy, float error_scale, bool hw,
        bool bin_vis, std::uint32_t num_clusters, std::uint32_t samples_per_slice, float max_sample_perc, float ver_perc, bool show_vmat, bool ge);
    MatrixReconstructionRenderer(const MatrixReconstructionRenderer& other) = delete;
    MatrixReconstructionRenderer(MatrixReconstructionRenderer&& other);
    MatrixReconstructionRenderer& operator = (const MatrixReconstructionRenderer& other) = delete;
    MatrixReconstructionRenderer& operator = (MatrixReconstructionRenderer&& other);
    ~MatrixReconstructionRenderer();

    bool render(Scene* scene, std::uint32_t spp, const RenderJob *job);

    void setCancel(bool cancel){
        std::lock_guard<std::mutex> lock(cancel_lock_);

        cancel_ = cancel;
    }

private:
    std::tuple<std::uint64_t, std::uint64_t> renderNonHW(Scene* scene, std::uint32_t spp, const RenderJob *job,
        std::vector<float>& timings, const std::vector<KDTNode<ReconstructionSample>*>& slices, 
        std::uint32_t samples_per_slice);
    std::tuple<std::uint64_t, std::uint64_t> renderHW(Scene* scene, std::uint32_t spp, const RenderJob *job,
        std::vector<float>& timings, const std::vector<KDTNode<ReconstructionSample>*>& slices, 
        std::uint32_t samples_per_slice, std::uint32_t slice_size);

private:
    std::vector<VPL> vpls_;
    float sample_percentage_, min_dist_, step_size_factor_, tolerance_, tau_;
    std::uint32_t max_iterations_, slice_size_;
    bool output_stats_;
    bool visibility_only_;
    bool adaptive_col_sampling_;
    bool adaptive_importance_sampling_;
    bool adaptive_force_resample_;
    bool adaptive_recover_transpose_;
    bool truncated_;
    bool show_slices_;
    bool vsl_;
    bool gather_stat_images_;
    bool show_svd_;
    ClusteringStrategy clustering_strategy_;
    float error_scale_;
    bool hw_;
    bool bin_vis_;
    std::uint32_t num_clusters_;
    std::uint32_t samples_per_slice_;
    float sample_inc_, max_sample_perc_;
    bool show_vmat_;
    bool ge_;
    std::mutex cancel_lock_;
    bool cancel_;

    std::vector<ReconstructionSample> samples_;
};

MTS_NAMESPACE_END

#endif