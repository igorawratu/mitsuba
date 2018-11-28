#ifndef PASSTHROUGHCLUSTERER_H_
#define PASSTHROUGHCLUSTERER_H_

#include <vector>

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

#include "manylightsbase.h"

MTS_NAMESPACE_BEGIN

class PassthroughClusterer : public ManyLightsClusterer{
public:
    PassthroughClusterer();
    PassthroughClusterer(const std::vector<VPL>& vpls) : vpls_(vpls){}
    PassthroughClusterer(const PassthroughClusterer& other) : vpls_(other.vpls_){}
    PassthroughClusterer(PassthroughClusterer&& other) : vpls_(std::move(other.vpls_)){}
    PassthroughClusterer& operator = (const PassthroughClusterer& other){
        if(this != &other){
            vpls_ = other.vpls_;
        }

        return *this;
    }

    PassthroughClusterer& operator = (PassthroughClusterer&& other){
        if(this != &other){
            vpls_ = std::move(other.vpls_);
        }

        return *this;
    }

    ~PassthroughClusterer(){}

    std::vector<VPL> getClusteringForPoint(const Intersection& its){
        return vpls_;
    }

private:
    std::vector<VPL> vpls_;
};

MTS_NAMESPACE_END

#endif