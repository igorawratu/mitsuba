#ifndef LIGHTTREE_H_
#define LIGHTTREE_H_

#include <memory>

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>


struct LightTreeNode{
    LightTreeNode() : left(nullptr), right(nullptr), vpl(nullptr), emission_scale(0.f){
    }

    LightTreeNode(const LightTreeNode& other) : left(other.left == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.left)),
        right(other.right == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.right)),
        vpl(other.vpl), emission_scale(other.emission_scale), idx(other.idx){
    }

    LightTreeNode(LightTreeNode&& other) : left(std::move(other.left)), right(std::move(other.right)),
        vpl(other.vpl), emission_scale(other.emission_scale), idx(other.idx){
    }

    LightTreeNode& operator = (const LightTreeNode& other){
       if(this != &other){
           left = other.left == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.left);
           right = other.right == nullptr ? nullptr : std::make_unique<LightTreeNode>(other.right);
           vpl = other.vpl;
           emission_scale = other.emission_scale;
           idx = other.idx;
       }

       return *this;
    }

    LightTreeNode& operator = (LightTreeNode&& other){
        if(this != &other){
           left = other.left == std::move(other.left);
           right = other.right == std::move(other.right);
           vpl = other.vpl;
           emission_scale = other.emission_scale;
           idx = other.idx;
       }

       return *this;
    }

    std::unique_ptr<LightTreeNode> left, right;
    VPL* vpl;
    float emission_scale;
    std::uint32_t idx;
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
    float *point_similarity_matrix_, *directional_similarity_matrix_, *oriented_similarity_matrix_;
    std::unique_ptr<LightTreeNode> point_tree_root_, directional_tree_root_, oriented_tree_root_;
}

#endif