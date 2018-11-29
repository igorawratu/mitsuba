#include "matrixseparationrenderer.h"

#include <vector>
#include <limits>
#include <utility>

MTS_NAMESPACE_BEGIN

struct RowSample{
    RowSample(Vector3f p, Vector3f n, std::uint32_t x, std::uint32_t y, std::uint32_t cols) : position(p), normal(n),
        image_x(x), image_y(y), col_samples(cols){
    }

    Vector3f position, normal;
    std::uint32_t image_x, image_y;
    std::vector<Spectrum> col_samples;
};

struct KDTNode{
    KDTNode() : samples(), left(nullptr), right(nullptr){
    }

    std::vector<RowSample> samples;
    std::unique_ptr<KDTNode> left;
    std::unique_ptr<KDTNode> right;

    void Split(float norm_scale){
        float maxf = std::numeric_limits<float>::max();
        float minf = std::numeric_limits<float>::min();
        Vector3f min_pos(maxf, maxf, maxf), max_pos(minf, minf, minf);
        Vector3f min_normal(maxf, maxf, maxf), max_normal(minf, minf, minf);

        for(size_t i = 0; i < samples.size(); ++i){
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            min_pos.x = std::min(samples[i].position.x, min_pos.x);
            max_pos.x = std::min(samples[i].position.x, max_pos.x);
            max_pos.x = std::min(samples[i].position.x, max_pos.x);
            max_pos.x = std::min(samples[i].position.x, max_pos.x);

            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            min_normal.x = std::min(samples[i].normal.x, min_normal.x);
            max_normal.x = std::min(samples[i].normal.x, max_normal.x);
            max_normal.x = std::min(samples[i].normal.x, max_normal.x);
            max_normal.x = std::min(samples[i].normal.x, max_normal.x);
        }

        std::array<std::pair<std::uint8_t, float>, 6> ranges;
        ranges[0] = std::make_pair(0, max_pos.x - min_pos.x);
        ranges[1] = std::make_pair(1, max_pos.y - min_pos.y);
        ranges[2] = std::make_pair(2, max_pos.z - min_pos.z);
        ranges[3] = std::make_pair(3, (max_normal.x - min_normal.x) * norm_scale);
        ranges[4] = std::make_pair(4, (max_normal.x - min_normal.x) * norm_scale);
        ranges[5] = std::make_pair(5, (max_normal.x - min_normal.x) * norm_scale);

        std::sort(ranges.begin(), ranges.end(), 
            [](const std::pair<std::uint8_t, float>& lhs, const std::pair<std::uint8_t, float>& rhs){
                return lhs.second > rhs.second;
            });

        switch(ranges[0].first){
            case 0:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.x > rhs.position.x;
                });
                break;
            case 1:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.y > rhs.position.y;
                });
                break;
            case 2:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.z > rhs.position.z;
                });
                break;
            case 3:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.normal.x > rhs.normal.x;
                });
                break;
            case 4:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.normal.y > rhs.normal.y;
                });
                break;
            case 5:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.normal.z > rhs.normal.z;
                });
                break;
            default:
                std::sort(samples.begin(), samples.end(), [](const RowSample& lhs, const RowSample& rhs){
                    return lhs.position.x > rhs.position.x;
                });
                break;
        }

        left = std::make_unique<KDTNode>();
        right = std::make_unique<KDTNode>();

        std::copy(samples.begin(), samples.begin() + samples.size() / 2, left->samples.begin());
        std::copy(samples.begin() + samples.size() / 2, samples.end(), right->samples.begin());
    }
};



bool MatrixSeparationRenderer::Render(Scene* scene){

}

MTS_NAMESPACE_END