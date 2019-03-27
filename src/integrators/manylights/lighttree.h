#ifndef LIGHTTREE_H_
#define LIGHTTREE_H_

#include <memory>
#include <mutex>

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

struct LightTreeNode{
    LightTreeNode() : left(nullptr), right(nullptr), vpl(EPointEmitterVPL, Spectrum(0.f)), emission_scale(0.f), 
        min_bounds(0.f), max_bounds(0.f), cone_ray(0.f), cone_halfangle(0.f){
    }

    LightTreeNode(const LightTreeNode& other) : left(other.left == nullptr ? nullptr : new LightTreeNode(*other.left)),
        right(other.right == nullptr ? nullptr : new LightTreeNode(*other.right)),
        vpl(other.vpl), emission_scale(other.emission_scale), min_bounds(other.min_bounds), max_bounds(other.max_bounds), 
        cone_ray(other.cone_ray), cone_halfangle(other.cone_halfangle){
    }

    LightTreeNode(LightTreeNode&& other) : left(std::move(other.left)), right(std::move(other.right)),
        vpl(other.vpl), emission_scale(other.emission_scale), min_bounds(other.min_bounds), max_bounds(other.max_bounds), 
        cone_ray(other.cone_ray), cone_halfangle(other.cone_halfangle){
    }

    LightTreeNode& operator = (const LightTreeNode& other){
       if(this != &other){
           left = other.left == nullptr ? nullptr : std::unique_ptr<LightTreeNode>(new LightTreeNode(*other.left));
           right = other.right == nullptr ? nullptr : std::unique_ptr<LightTreeNode>(new LightTreeNode(*other.right));
           vpl = other.vpl;
           emission_scale = other.emission_scale;
           min_bounds = other.min_bounds;
           max_bounds = other.max_bounds;
           cone_ray = other.cone_ray;
           cone_halfangle = other.cone_halfangle;
       }

       return *this;
    }

    LightTreeNode& operator = (LightTreeNode&& other){
        if(this != &other){
           left = std::move(other.left);
           right = std::move(other.right);
           vpl = other.vpl;
           emission_scale = other.emission_scale;
           min_bounds = other.min_bounds;
           max_bounds = other.max_bounds;
           cone_ray = other.cone_ray;
           cone_halfangle = other.cone_halfangle;

           other.left = nullptr;
           other.right = nullptr;
       }

       return *this;
    }

    std::unique_ptr<LightTreeNode> left, right;
    VPL vpl;
    float emission_scale;
    Point min_bounds, max_bounds;
    Vector3 cone_ray;
    float cone_halfangle;
};

class LightTree : public ManyLightsClusterer{
public:
    LightTree();
    LightTree(const std::vector<VPL>& vpls, float min_dist, std::uint32_t max_lights, float error_threshold);
    LightTree(const LightTree& other);
    LightTree(LightTree&& other);
    LightTree& operator = (const LightTree& other);
    LightTree& operator = (LightTree&& other);
    ~LightTree();

    std::vector<VPL> getClusteringForPoint(const Intersection& its);
    std::vector<VPL> getClusteringForPoints(const std::vector<Intersection>& points);

private:
    std::vector<VPL> point_vpls_, directional_vpls_, oriented_vpls_;
    std::unique_ptr<LightTreeNode> point_tree_root_, directional_tree_root_, oriented_tree_root_;
    float min_dist_;
    std::uint32_t max_lights_;
    float error_threshold_;
    std::mutex mutex_;
};

MTS_NAMESPACE_END

#endif