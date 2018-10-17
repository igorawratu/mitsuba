#ifndef MANYLIGHTSRENDERER_H_
#define MANYLIGHTSRENDERER_H_

template <class Clusterer, class Renderer>
class ManyLightsRenderer{
public:
    ManyLightsRenderer();
    ManyLightsRenderer(const ManyLightsRenderer& other);
    ManyLightsRenderer(ManyLightsRenderer&& other);
    ManyLightsRenderer& operator = (const ManyLightsRenderer& other);
    ManyLightsRenderer& operator = (ManyLightsRenderer&& other);
    ~ManyLightsRenderer();

    void Render();

private:
    Clusterer clusterer_;
    Renderer renderer_;
};

#endif