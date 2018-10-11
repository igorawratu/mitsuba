#include "rowcolumnsampling.h"

#include <random>
#include <queue>
#include <pair>

//works by assuming a z-curve ordering so to reduce the bucketing and samplingproblem to 1D, then reverses the 
//z-ordering process to acquire the x and y values, and by extension the actual index of the pixel which has a 
//row-col ordering
std::vector<std::pair<std::uint32_t, std::uint32_t>> subsampleRows(std::uint32_t rows, std::tuple<std::uint32_t, std::uint32_t> resolution){
    std::vector<std::pair<std::uint32_t, std::uint32_t>> row_indices;

    //find factors to try get a decent bucket size for the requested number of buckets.   
    typedef std::pair<std::int32_t, std::int32_t> FactorPair;

    auto comparator = [](FactorPair l, FactorPair r){
        return abs(l.first - l.second) > abs(r.first - r.second);
    };

    std::priority_queue<FactorPair, std::vector<FactorPair>, decltype(comparator)> factors(comparator);

    for(int i = 0; i < (int)(sqrt((float)rows) + 1.f); ++i){
        if(rows % i == 0){
            factors.push_back(std::make_pair(i, rows / i));
        }
    }

    auto dims = factors.top();
    bool rows_more = std::get<0>(resolution) > std::get<1>(resolution);
    std::uint32_t rows = rows_more ? std::max(dims.first, dims.second) : std::min(dims.first, dims.second);
    std::uint32_t cols = rows_more ? std::min(dims.first, dims.second) : std::max(dims.first, dims.second);

    std::uint32_t bucket_width = std::get<0>(resolution) / rows;
    std::uint32_t bucket_height = std::get<1>(resolution) / cols;

    rows += std::get<0>(resolution) % rows == 0 ? 0 : 1;
    cols += std::get<1>(resolution) % cols == 0 ? 0 : 1;

    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);
    
    for(std::uint32_t i = 0; i < rows; ++i){
        for(std::uint32_t j = 0; j < cols; ++j){
            std::uint32_t x_start = j * bucket_width;
            std::uint32_t y_start = i * bucket_height;

            auto gen_y = std::bind(std::uniform_int_distribution<std::uint32_t>(0, 
                        std::min(bucket_height, std::get<1>(resolution) - y_start) - 1));
            auto gen_x = std::bind(std::uniform_int_distribution<std::uint32_t>(0, 
                        std::min(bucket_width, std::get<0>(resolution) - x_start) - 1));

            int x = gen_x() + x_start;
            int y = gen_y() + y_start;

            row_indices.push_back(std::make_pair(x, y));
        }
    }

    return row_indices;
}

std::vector<float> calculateLightContributions(const std::vector<VPL>& vpls, 
                                    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& rows,
                                    Scene* scene, float min_dist){
    std::vector<float> contributions(vpls.size());

    ref<Sensor> sensor = scene->getSensor();

    //disregarding aperture and time sampling for now, as we are only dealing with a single sample per pixel
    Point2 aperture_sample(0.5f, 0.5f);
    Float time_sample(0.5f);

    for(int i = 0; i < vpls.size(); ++i){
        //only dealing with emitter and surface VPLs curently.
        if (vpls[i].type != EPointEmitterVPL && vpls[i].type != ESurfaceVPL){
            continue;
        }

        for(int j = 0; j < rows.size(); ++j){
            Point2 sample_position(rows[j].first + 0.5f, rows[j].second + 0.5f);

            Ray ray;
            sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

            Intersection its;

            if (!scene->rayIntersect(ray, its)) {
                continue;
            }

            Normal n = its.geoFrame.n;

            Point ray_origin = its.p;
            Ray shadow_ray(ray_origin, normalize(vpls[i].its.p - ray_origin), ray.time);

            Float t;
            ConstShapePtr shape;
            Normal norm;
            Point2 uv;

            if(scene->rayIntersect(shadow_ray, t, shape, norm, uv)){
                if(abs((ray_origin - vpls[i].its.p).length() - t) > 0.0001f ){
                    continue;
                }
            }

            BSDFSamplingRecord bsdf_sample_record(its, sampler, ERadiance);
            bsdf_sample_record.wi = its.toLocal(normalize(vpls[i].its.p - its.p));
            bsdf_sample_record.wo = its.toLocal(n);

            albedo = its.getBSDF()->eval(bsdf_sample_record);

            float d = std::max((its.p - vpls[i].its.p).length(), min_dist);
            float attenuation = 1.f / (d * d);

            float n_dot_ldir = std::max(0.f, dot(normalize(n), normalize(vpls[i].its.p - its.p)));
            float ln_dot_ldir = std::max(0.f, dot(normalize(vpls[i].its.shFrame.n), normalize(its.p - vpls[i].its.p)));

            Spectrum lightContribution = (vpls[i].P * ln_dot_ldir * attenuation * n_dot_ldir * albedo) / PI;
            contributions[i] += lightContribution.getLuminance();
        }
    }

    return contributions;
}

std::vector<VPL> calculateClustering(const std::vector<VPL>& vpls, const std::vector<float>& contributions, std::uint32_t num_clusters){
    std::uint32_t clusters_by_sampling = num_clusters > 1 ? num_clusters * 2 / 3 : num_clusters;

    std::vector<std::vector<VPL>> clusters;

    while(clusters.size() < clusters_by_sampling){
        
    }
}

RowColumnSampling::RowColumnSampling(const std::vector<VPL>& vpls, std::uint32_t num_clusters, std::uint32_t rows, std::uint32_t cols,
        std::tuple<std::uint32_t, std::uint32_t> resolution, Scene* scene, float min_dist){


}

RowColumnSampling::RowColumnSampling(const RowColumnSampling& other){

}

RowColumnSampling::RowColumnSampling(RowColumnSampling&& other){

}

RowColumnSampling::RowColumnSampling& operator = (const RowColumnSampling& other){

}

RowColumnSampling::RowColumnSampling& operator = (RowColumnSampling&& other){

}

RowColumnSampling::~RowColumnSampling(){

}

void RowColumnSampling::updateClustering(const std::vector<VPL>& vpls, std::uint32_t num_clusters, std::uint32_t rows,
        std::uint32_t cols, std::tuple<std::uint32_t, std::uint32_t> resolution, Scene* scene, float min_dist){

}

std::vector<VPL> RowColumnSampling::getClusteringForPoint(){

}