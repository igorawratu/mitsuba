#ifndef MANYLIGHTSRENDERER_H_
#define MANYLIGHTSRENDERER_H_

#include <tuple>

#include <mitsuba/render/vpl.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

//needs this otherwise the static_assert evaluates immediately rather than when the incorrect templated class is
//actually instantiated
template<typename... Ts> struct fake_dependency: public std::false_type {};

template<typename... Ts>
class ManyLightsRenderer{
    ManyLightsRenderer(){
        static_assert(fake_dependency<Ts...>::value, "Many lights renderer specialization must be in the form of <ClustererType, std::tuple<ClusterInitParams...>, RendererType, std::tuple<RendererInitParams...>>");
    }
};

template<typename ClustererType, typename... ClusterInitParams, typename RendererType, typename... RendererInitParams> 
class ManyLightsRenderer<ClustererType, std::tuple<ClusterInitParams...>, RendererType, std::tuple<RendererInitParams...>>{
public:  
    ManyLightsRenderer(ClusterInitParams... cluster_init_params, RendererInitParams... renderer_init_params) : 
        clusterer_(cluster_init_params...), renderer_(renderer_init_params...){
    }

    void render(Scene* scene){
        renderer_.render(clusterer_, scene);
    }

    void setCancel(bool cancel_){
        renderer_.setCancel(cancel_);
    }

private:
    ClustererType clusterer_;
    RendererType renderer_;
};

MTS_NAMESPACE_END

#endif