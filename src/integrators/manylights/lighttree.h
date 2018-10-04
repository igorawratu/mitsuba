#ifndef LIGHTTREE_H_
#define LIGHTTREE_H_

#include <memory>

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

struct LightTreeNode{
    LightTreeNode() : left(nullptr), right(nullptr), emission_scale(0.f), tl(0.f), tr(0.f), cone_ray1(0.f),
        cone_ray2(0.f){
    }

    LightTreeNode(const LightTreeNode& other) : left(other.left == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.left)),
        right(other.right == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.right)),
        vpl(other.vpl), emission_scale(other.emission_scale), tl(other.tl), br(other.br), cone_ray1(other.cone_ray1),
        cone_ray2(other.cone_ray2){
    }

    LightTreeNode(LightTreeNode&& other) : left(std::move(other.left)), right(std::move(other.right)),
        vpl(other.vpl), emission_scale(other.emission_scale), tl(other.tl), br(other.br), cone_ray1(other.cone_ray1),
        cone_ray2(other.cone_ray2){
    }

    LightTreeNode& operator = (const LightTreeNode& other){
       if(this != &other){
           left = other.left == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.left);
           right = other.right == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.right);
           vpl = other.vpl;
           emission_scale = other.emission_scale;
           tl = other.tl;
           br = other.br;
           cone_ray1 = other.cone_ray1;
           cone_ray2 = other.cone_ray2;
       }

       return *this;
    }

    LightTreeNode& operator = (LightTreeNode&& other){
        if(this != &other){
           left = std::move(other.left);
           right = std::move(other.right);
           vpl = other.vpl;
           emission_scale = other.emission_scale;
           tl = other.tl;
           br = other.br;
           cone_ray1 = other.cone_ray1;
           cone_ray2 = other.cone_ray2;
       }

       return *this;
    }

    std::unique_ptr<LightTreeNode> left, right;
    VPL vpl;
    float emission_scale;
    Point tl, br;
    Vector3 cone_ray1, cone_ray2;
};

class LightTree{
public:
    LightTree();
    LightTree(const std::vector<VPL>& vpls);
    LightTree(const LightTree& other);
    LightTree(const LightTree&& other);
    LightTree& operator = (const LightTree& other);
    LightTree& operator = (LightTree&& other);
    ~LightTree();

    void setVPLs(const std::vector<VPL>& vpls);
    std::vector<VPL> getClusteringForPoint(Intersection its);

private:
    std::vector<VPL> point_vpls_, directional_vpls_, oriented_vpls_;
    std::unique_ptr<LightTreeNode> point_tree_root_, directional_tree_root_, oriented_tree_root_;
};

MTS_NAMESPACE_END

#endif