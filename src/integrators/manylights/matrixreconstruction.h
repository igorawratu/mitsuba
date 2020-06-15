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
        image_x(x), image_y(y), intersected_scene(in_scene), its(intersection), ray(r), color(0.f){
    }

    ReconstructionSample(ReconstructionSample&& other) : image_x(other.image_x), image_y(other.image_y), 
        intersected_scene(other.intersected_scene), its(other.its), ray(other.ray), color(other.color),
        unoccluded_samples(other.unoccluded_samples){
    }

    ReconstructionSample& operator = (ReconstructionSample&& other){
        if(this != &other){
            image_x = other.image_x;
            image_y = other.image_y;
            intersected_scene = other.intersected_scene;
            its = other.its;
            ray = other.ray;
            color = other.color;
            unoccluded_samples = other.unoccluded_samples;
        }

        return *this;
    }

    std::uint32_t image_x, image_y;
    bool intersected_scene;
    Intersection its;
    Ray ray;
    Spectrum color;
    std::vector<Spectrum> unoccluded_samples;
};

enum ClusteringStrategy{MDLC, LS};

class MatrixReconstructionRenderer : public ManyLightsRenderer{
public:
    MatrixReconstructionRenderer() = delete;
    MatrixReconstructionRenderer(const std::vector<VPL>& vpls, float sample_percentage, 
        float min_dist, std::uint32_t slice_size, bool adaptive_importance_sampling, 
        bool show_slices, bool vsl, ClusteringStrategy clustering_strategy, bool hw,
        bool bin_vis, std::uint32_t num_clusters, std::uint32_t samples_per_slice, float max_sample_perc, float ver_perc, bool ge);
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
    float sample_percentage_, min_dist_;
    std::uint32_t slice_size_;
    bool adaptive_importance_sampling_;
    bool show_slices_;
    bool vsl_;
    ClusteringStrategy clustering_strategy_;
    bool hw_;
    bool bin_vis_;
    std::uint32_t num_clusters_;
    std::uint32_t samples_per_slice_;
    float sample_inc_, max_sample_perc_;
    bool ge_;
    std::mutex cancel_lock_;
    bool cancel_;

    std::vector<ReconstructionSample> samples_;
};

MTS_NAMESPACE_END

#endif