#ifndef ILLUMCUT_H_
#define ILLUMCUT_H_

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/statistics.h>
#include <vector>
#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

#include "lighttree.h"

MTS_NAMESPACE_BEGIN

struct IllumcutSample{
    IllumcutSample() : intersected_scene(false){
    }

    IllumcutSample(std::uint32_t x, std::uint32_t y, bool in_scene, Intersection intersection, const Ray& r) : 
        image_x(x), image_y(y), intersected_scene(in_scene), its(intersection), ray(r), color(0.f), upper_bound(0.f){
    }

    IllumcutSample(IllumcutSample&& other) : image_x(other.image_x), image_y(other.image_y), 
        intersected_scene(other.intersected_scene), its(other.its), ray(other.ray), color(other.color), upper_bound(other.upper_bound){
    }

    IllumcutSample& operator = (IllumcutSample&& other){
        if(this != &other){
            image_x = other.image_x;
            image_y = other.image_y;
            intersected_scene = other.intersected_scene;
            its = other.its;
            ray = other.ray;
            color = other.color;
            upper_bound = other.upper_bound;
        }

        return *this;
    }

    std::uint32_t image_x, image_y;
    bool intersected_scene;
    Intersection its;
    Ray ray;
    Spectrum color;

    float upper_bound;
};

class IlluminationCutRenderer : public ManyLightsRenderer{
public:
    IlluminationCutRenderer(const std::vector<VPL>& vpls, float error_thresh, float min_dist, float upper_distance_thresh,
        std::uint32_t num_clusters);
    IlluminationCutRenderer(const IlluminationCutRenderer& other) = delete;
    IlluminationCutRenderer(IlluminationCutRenderer&& other) = delete;
    IlluminationCutRenderer& operator = (const IlluminationCutRenderer& other) = delete;
    IlluminationCutRenderer& operator = (IlluminationCutRenderer&& other) = delete;
    ~IlluminationCutRenderer();

    bool render(Scene* scene, std::uint32_t spp, const RenderJob *job);
    void setCancel(bool cancel){}

private:
    std::vector<VPL> vpls_;
    float error_threshold_, min_dist_;
    std::vector<IllumcutSample> samples_;
    float upper_distance_thresh_;
    std::uint32_t num_clusters_;
};

MTS_NAMESPACE_END

#endif