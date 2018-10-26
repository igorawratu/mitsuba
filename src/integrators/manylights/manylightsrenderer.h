#ifndef MANYLIGHTSRENDERER_H_
#define MANYLIGHTSRENDERER_H_

template <class Clusterer, class Renderer>
class ManyLightsRenderer{
public:
    ManyLightsRenderer() : clusterer_(), renderer_(){}
    ManyLightsRenderer(const ManyLightsRenderer& other) : clusterer_(other.clusterer_), renderer_(other.renderer_){}
    ManyLightsRenderer(ManyLightsRenderer&& other) : clusterer_(std::move(other.clusterer_)), renderer_(std::move(other.renderer_)){}
    ManyLightsRenderer& operator = (const ManyLightsRenderer& other);
    ManyLightsRenderer& operator = (ManyLightsRenderer&& other);
    ~ManyLightsRenderer();

    void Render();
    void BindOutputbuffer(void* buffer);

private:
    Clusterer clusterer_;
    Renderer renderer_;
};

#endif