#include "illumcut.h"

#include <mutex>
#include <memory>
#include <queue>
#include <tuple>
#include <stack>
#include <chrono>
#include <random>
#include <unordered_map>
#include <iostream>
#include <deque>
#include <vector>

#include "common.h"
#include "filewriter.h"

MTS_NAMESPACE_BEGIN

IlluminationCutRenderer::IlluminationCutRenderer(const std::vector<VPL>& vpls, float error_thresh, float min_dist, float upper_distance_thresh,
    std::uint32_t num_clusters) :
    vpls_(vpls),
    error_threshold_(error_thresh),
    min_dist_(min_dist),
    upper_distance_thresh_(upper_distance_thresh),
    num_clusters_(num_clusters){

}


IlluminationCutRenderer::~IlluminationCutRenderer(){

}


void expandBB(std::pair<Vector3f, Vector3f>& bb, const Vector3f& v){
    bb.first.x = std::min(bb.first.x, v.x);
    bb.first.y = std::min(bb.first.y, v.y);
    bb.first.z = std::min(bb.first.z, v.z);

    bb.second.x = std::max(bb.second.x, v.x);
    bb.second.y = std::max(bb.second.y, v.y);
    bb.second.z = std::max(bb.second.z, v.z);
}

std::unique_ptr<OctreeNode<IllumcutSample>> constructOctree(Scene* scene, std::vector<IllumcutSample>& samples, float min_dist, std::uint32_t spp){
    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    samples.resize(film->getSize().y * film->getSize().x * spp);

    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    for (std::int32_t y = 0; y < film->getSize().y; ++y) {
        for (std::int32_t x = 0; x < film->getSize().x; ++x) {
            std::uint32_t cell_dim = sqrt(spp) + 0.5f;
            float cell_side_len = 1.f / cell_dim;

            for(std::uint32_t i = 0; i < spp; ++i){
                //calculates the position the ray intersects with
                Ray ray;

                float x_jitter = 0.5f * cell_side_len;
                float y_jitter = 0.5f * cell_side_len;
                float x_off = (i % cell_dim) * cell_side_len;
                float y_off = (i / cell_dim) * cell_side_len;

                Point2 sample_position(x + x_off + x_jitter, y + y_off + y_jitter);

                Point2 aperture_sample(0.5f, 0.5f);
                Float time_sample(0.5f);

                sensor->sampleRay(ray, sample_position, aperture_sample, time_sample);

                IllumcutSample curr_sample;
                curr_sample.image_x = x;
                curr_sample.image_y = y;
                curr_sample.ray = ray;
                curr_sample.color = Spectrum(0.f);

                std::uint32_t num_bounces = 0;

                while(true){
                    curr_sample.intersected_scene = scene->rayIntersect(curr_sample.ray, curr_sample.its);
                    if(!curr_sample.intersected_scene){
                        break;
                    }

                    if((curr_sample.its.getBSDF()->getType() & BSDF::ESmooth) || curr_sample.its.isEmitter()
                        || ++num_bounces > 10){
                        break;
                    }

                    BSDFSamplingRecord bsdf_sample_record(curr_sample.its, sampler);
                    bsdf_sample_record.typeMask = curr_sample.its.getBSDF()->isDielectric() ? 
                        BSDF::EDeltaTransmission | BSDF::ENull : BSDF::EDeltaReflection;
                    
                    curr_sample.its.getBSDF()->sample(bsdf_sample_record, sampler->next2D());

                    curr_sample.ray = Ray(curr_sample.its.p, bsdf_sample_record.its.toWorld(bsdf_sample_record.wo), ray.time);
                }

                if(!curr_sample.intersected_scene){
                    if(scene->hasEnvironmentEmitter()){
                        curr_sample.color = scene->evalEnvironment(RayDifferential(curr_sample.ray));
                    }
                    else curr_sample.color = Spectrum(0.f);
                }
                else if(curr_sample.its.isEmitter()){
                    curr_sample.color = curr_sample.its.Le(-curr_sample.ray.d);
                }

                samples[y * film->getSize().x * spp + x * spp + i] = std::move(curr_sample);
            }
        }
    }

    //only samples that intersected the scene are recovered. The rest are sampled from the environment map if there is one
    std::vector<std::uint32_t> sample_indices;
    sample_indices.reserve(samples.size());
    
    std::pair<Vector3f, Vector3f> bb(Vector3f(std::numeric_limits<float>::max()), Vector3f(-std::numeric_limits<float>::max()));
    std::pair<Vector3f, Vector3f> nbb(Vector3f(std::numeric_limits<float>::max()), Vector3f(-std::numeric_limits<float>::max()));

    for(std::uint32_t i = 0; i < samples.size(); ++i){
        if(samples[i].intersected_scene){
            sample_indices.push_back(i);
            Vector3f n(samples[i].its.shFrame.n);
            Vector3f p(samples[i].its.p);

            expandBB(bb, p);
            expandBB(nbb, n);
        }
    }

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    return std::unique_ptr<OctreeNode<IllumcutSample>>(new OctreeNode<IllumcutSample>(&samples, std::move(sample_indices), bb, nbb, 0, 3, rng));
}

typedef std::pair<LightTreeNode*, OctreeNode<IllumcutSample>*> IllumPair;
const float CONE_THRESH = cos(90.f / 180.f * M_PI);

std::mutex printmut;

bool refineUpper(const IllumPair& illum_pair, float upper_distance_thresh){
    bool refine = false;

    float r1, r2, d;

    if(illum_pair.first->vpl.type != EDirectionalEmitterVPL){
        Vector3f c1 = Vector3f(illum_pair.first->min_bounds + illum_pair.first->max_bounds) / 2.f;
        Vector3f dim1 = Vector3f(illum_pair.first->max_bounds - illum_pair.first->min_bounds);
        r1 = std::max(dim1.x, std::max(dim1.y, dim1.z));

        Vector3f c2 = (illum_pair.second->bb.first + illum_pair.second->bb.second) / 2.f;
        Vector3f dim2 = illum_pair.second->bb.second - illum_pair.second->bb.first;
        r2 = std::max(dim2.x, std::max(dim2.y, dim2.z));

        d = std::max(0.f, upper_distance_thresh * ((c1 - c2).length()));

        refine |= std::max(r1, r2) > d;

    }

    if(illum_pair.first->vpl.type != EPointEmitterVPL){
        refine |= illum_pair.first->bcone.GetAngleCos() < CONE_THRESH;
    }

    return refine;   
}

bool isIllumAware(const IllumPair& illum_pair, float min_dist, float error_thresh){
    float max_error = 0.f;

    for(std::uint32_t i = 0; i < illum_pair.second->sample_indices.size(); ++i){
        IllumcutSample& rep_sample = illum_pair.second->representative();
        max_error = std::max(LightTree::calculateClusterBounds(rep_sample.its.p, rep_sample.its.shFrame.n, illum_pair.first, 
            illum_pair.first->vpl.type, min_dist), max_error);
    }

    return max_error < error_thresh * illum_pair.second->upper_bound;
}

bool lines_overlap(float x1, float x2, float y1, float y2){
    bool overlapping = false;

    overlapping |= x1 <= y1 && y1 <= x2;
    overlapping |= x1 <= y2 && y2 <= x2;
    overlapping |= y1 <= x1 && x1 <= y2;
    overlapping |= y1 <= x2 && x2 <= y2;

    return overlapping;
}

float distanceSqr(const std::pair<Vector3f, Vector3f>& a, const std::pair<Vector3f, Vector3f>& b)
{
    float dist_sqr = 0;
    
    if(b.second.x < a.first.x)
    {
        float d = b.second.x - a.first.x;
        dist_sqr += d * d;
    }
    else if(b.first.x > a.second.x)
    {
        float d = b.first.x - a.second.x;
        dist_sqr += d * d;
    }

    if(b.second.y < a.first.y)
    {
        float d = b.second.y - a.first.y;
        dist_sqr += d * d;
    }
    else if(b.first.y > a.second.y)
    {
        float d = b.first.y - a.second.y;
        dist_sqr += d * d;
    }

    if(b.second.z < a.first.z)
    {
        float d = b.second.z - a.first.z;
        dist_sqr += d * d;
    }
    else if(b.first.z > a.second.z)
    {
        float d = b.first.z - a.second.z;
        dist_sqr += d * d;
    }

    return dist_sqr;
}

bool refineLTree(const IllumPair& illum_pair){
    if(illum_pair.first->left == nullptr || illum_pair.first->right == nullptr){
        return false;
    }

    if(illum_pair.second->sample_indices.size() == 1){
        return true;
    }

    assert(illum_pair.first->left && illum_pair.first->right);

    std::pair<Vector3f, Vector3f> lightbb(Vector3f(illum_pair.first->min_bounds), Vector3f(illum_pair.first->max_bounds));

    Vector3f lbb_extents = lightbb.second - lightbb.first;
    Vector3f rbb_extents = illum_pair.second->bb.second - illum_pair.second->bb.first;

    if(illum_pair.first->vpl.type != EDirectionalEmitterVPL &&
        (lines_overlap(lightbb.first.x, lightbb.second.x, illum_pair.second->bb.first.x, illum_pair.second->bb.second.x) &&
        lines_overlap(lightbb.first.y, lightbb.second.y, illum_pair.second->bb.first.y, illum_pair.second->bb.second.y) &&
        lines_overlap(lightbb.first.z, lightbb.second.z, illum_pair.second->bb.first.z, illum_pair.second->bb.second.z))){

        return lbb_extents.lengthSquared() > rbb_extents.lengthSquared();
    }
    
    
    float dsqr = illum_pair.first->vpl.type != EDirectionalEmitterVPL ? distanceSqr(lightbb, illum_pair.second->bb) : 1.f;

    float lg = 1.f;
    if(illum_pair.first->vpl.type != EDirectionalEmitterVPL){
        std::pair<Vector3f, Vector3f> lightc1bb(Vector3f(illum_pair.first->left->min_bounds), Vector3f(illum_pair.first->left->max_bounds));
        std::pair<Vector3f, Vector3f> lightc2bb(Vector3f(illum_pair.first->right->min_bounds), Vector3f(illum_pair.first->right->max_bounds));

        lg = std::min(distanceSqr(lightc1bb, illum_pair.second->bb), distanceSqr(lightc2bb, illum_pair.second->bb)) / dsqr;
    }

    float lm = lbb_extents.length() / sqrt(dsqr);

    if(illum_pair.first->vpl.type != EPointEmitterVPL){
        lm *= (1.0f - illum_pair.first->bcone.GetAngleCos());
    }
    
    float lheuristic = lm * lg;

    float rg = 1.f;
    if(illum_pair.first->vpl.type != EDirectionalEmitterVPL){
        rg = std::numeric_limits<float>::max();

        for(std::uint32_t i = 0; i < illum_pair.second->children.size(); ++i){
            if(illum_pair.second->children[i]){
                rg = std::min(rg, distanceSqr(lightbb, illum_pair.second->children[i]->bb));
            }
        }

        rg /= dsqr;
    }

    float rm = (1.0f - illum_pair.second->bcone.GetAngleCos()) * rbb_extents.length() / sqrt(dsqr);

    float rheuristic = rm * rg;

    bool use_heuristics = std::abs(lheuristic - rheuristic) > 0.0001f;

    if(use_heuristics){
        return lheuristic > rheuristic;
    }
    else{
        return lbb_extents.lengthSquared() > rbb_extents.lengthSquared();
    }
      
}

typedef std::tuple<LightTreeNode*, float> ClusterAndScore;

void getAllLeaves(OctreeNode<IllumcutSample>* node, std::vector<OctreeNode<IllumcutSample>*>& leaves){
    if(node->sample_indices.size() == 1){
        leaves.push_back(node);
    }
    else{
        for(std::uint32_t i = 0; i < node->children.size(); ++i){
            if(node->children[i] != nullptr){
                getAllLeaves(node->children[i].get(), leaves);
            }
        }
    }
}

void computeUpperBounds2Helper(LightTree* lt, std::vector<OctreeNode<IllumcutSample>*> leaves, float min_dist, std::uint32_t num_clusters){
    auto comparator = [](ClusterAndScore l, ClusterAndScore r){
        return std::get<1>(l) < std::get<1>(r);
    };

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < leaves.size(); ++i){
        //largest first, but front of the queue is last out
        
        std::priority_queue<ClusterAndScore, std::vector<ClusterAndScore>, decltype(comparator)> pqueue(comparator);

        IllumcutSample& curr_sample = leaves[i]->representative();

        if(lt->getPointTreeRoot() != nullptr){
            float error = LightTree::calculateClusterBounds(curr_sample.its.p, curr_sample.its.shFrame.n, lt->getPointTreeRoot(), 
                    lt->getPointTreeRoot()->vpl.type, min_dist);
            pqueue.push(std::make_tuple(lt->getPointTreeRoot(), error));
        }

        if(lt->getOrientedTreeRoot() != nullptr){
            float error = LightTree::calculateClusterBounds(curr_sample.its.p, curr_sample.its.shFrame.n, lt->getOrientedTreeRoot(), 
                    lt->getOrientedTreeRoot()->vpl.type, min_dist);
            pqueue.push(std::make_tuple(lt->getOrientedTreeRoot(), error));
        }

        if(lt->getDirectionalTreeRoot() != nullptr){
            float error = LightTree::calculateClusterBounds(curr_sample.its.p, curr_sample.its.shFrame.n, lt->getDirectionalTreeRoot(), 
                    lt->getDirectionalTreeRoot()->vpl.type, min_dist);
            pqueue.push(std::make_tuple(lt->getDirectionalTreeRoot(), error));
        }

        std::uint32_t num_added = 0;
        float bound = 0.f;

        while(pqueue.size() + num_added < num_clusters){
            ClusterAndScore entry = pqueue.top();
            pqueue.pop();

            LightTreeNode* node = std::get<0>(entry);
            
            if(node->left == nullptr && node->right == nullptr){
                num_added++;
                bound += std::get<1>(entry);
                continue;
            }

            if(node->left.get() != nullptr){
                float err = LightTree::calculateClusterBounds(curr_sample.its.p, curr_sample.its.shFrame.n, node->left.get(), node->left->vpl.type, min_dist);
                pqueue.push(std::make_tuple(node->left.get(), err));
            }

            if(node->right.get() != nullptr){
                float err = LightTree::calculateClusterBounds(curr_sample.its.p, curr_sample.its.shFrame.n, node->right.get(), node->right->vpl.type, min_dist);
                pqueue.push(std::make_tuple(node->right.get(), err));
            }
        }

        while(pqueue.size() > 0){
            auto entry = pqueue.top();
            pqueue.pop();

            bound += std::get<1>(entry);
        }

        leaves[i]->upper_bound = bound;
    }
}

void computeUpperBounds2(LightTree* lt, OctreeNode<IllumcutSample>* rt_root, float min_dist, std::uint32_t num_clusters){
    std::vector<OctreeNode<IllumcutSample>*> leaves;
    getAllLeaves(rt_root, leaves);

    computeUpperBounds2Helper(lt, leaves, min_dist, num_clusters);

    rt_root->cacheMinUpper();
}

void computeUpperBounds(LightTree* lt, OctreeNode<IllumcutSample>* rt_root, Scene* scene, float min_dist, float upper_distance_thresh){
    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    std::deque<IllumPair> initial_pairs;

    if(lt->getPointTreeRoot() != nullptr){
        initial_pairs.push_back(std::make_pair(lt->getPointTreeRoot(), rt_root));
    }

    if(lt->getDirectionalTreeRoot() != nullptr){
        initial_pairs.push_back(std::make_pair(lt->getDirectionalTreeRoot(), rt_root));
    }

    if(lt->getOrientedTreeRoot() != nullptr){
        initial_pairs.push_back(std::make_pair(lt->getOrientedTreeRoot(), rt_root));
    }

    while(initial_pairs.size() < 64){
        IllumPair curr = initial_pairs.front();
        initial_pairs.pop_front();

        for(std::uint8_t i = 0; i < curr.second->children.size(); ++i){
            if(curr.second->children[i] != nullptr){
                initial_pairs.push_back(std::make_pair(curr.first, curr.second->children[i].get()));
            }
        }
    }

    std::vector<IllumPair> initial_pairs_vec(initial_pairs.begin(), initial_pairs.end());

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < initial_pairs_vec.size(); ++i){
        std::stack<IllumPair> node_stack;
        node_stack.push(initial_pairs_vec[i]);

        while(!node_stack.empty()){
            IllumPair curr = node_stack.top();
            node_stack.pop();

            if(refineUpper(curr, upper_distance_thresh)){
                if(refineLTree(curr)){
                    if(curr.first->left != nullptr){
                        node_stack.push(std::make_pair(curr.first->left.get(), curr.second));
                    }

                    if(curr.first->right != nullptr){
                        node_stack.push(std::make_pair(curr.first->right.get(), curr.second));
                    }
                }
                else{
                    for(std::uint8_t i = 0; i < curr.second->children.size(); ++i){
                        if(curr.second->children[i] != nullptr){
                            node_stack.push(std::make_pair(curr.first, curr.second->children[i].get()));
                        }
                    }
                }
            }
            else{
                IllumcutSample& curr_sample = curr.second->representative();

                float estimated_error = LightTree::calculateClusterBounds(curr_sample.its.p, curr_sample.its.shFrame.n, curr.first, 
                    curr.first->vpl.type, min_dist);

                /*if(estimated_error < 0.00000001f){
                    std::lock_guard<std::mutex> lock(printmut);
                    Vector3f axis = curr.first->bcone.GetAxis();
                    std::cout << curr_sample.its.shFrame.n.x << " " << curr_sample.its.shFrame.n.y << " " << curr_sample.its.shFrame.n.z << "-" <<
                        axis.x << " " << axis.y << " " << axis.z << std::endl;
                }*/

                curr.second->updateUpperBound(estimated_error);
            }
        }
    }
    

    rt_root->cacheMinUpper();
}

std::vector<IllumPair> getIlluminationAwarePairs(LightTree* lt, OctreeNode<IllumcutSample>* rt_root, float min_dist, float error_thresh){
    
    std::vector<IllumPair> illum_aware_pairs;

    std::deque<IllumPair> initial_pairs;

    if(lt->getPointTreeRoot() != nullptr){
        initial_pairs.push_back(std::make_pair(lt->getPointTreeRoot(), rt_root));
    }

    if(lt->getDirectionalTreeRoot() != nullptr){
        initial_pairs.push_back(std::make_pair(lt->getDirectionalTreeRoot(), rt_root));
    }

    if(lt->getOrientedTreeRoot() != nullptr){
        initial_pairs.push_back(std::make_pair(lt->getOrientedTreeRoot(), rt_root));
    }

    while(initial_pairs.size() < 64){
        IllumPair curr = initial_pairs.front();
        initial_pairs.pop_front();
        
        for(std::uint8_t i = 0; i < curr.second->children.size(); ++i){
            if(curr.second->children[i] != nullptr){
                initial_pairs.push_back(std::make_pair(curr.first, curr.second->children[i].get()));
            }
        }
    }

    std::vector<IllumPair> initial_pairs_vec(initial_pairs.begin(), initial_pairs.end());

    std::mutex illum_aware_mutex;

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < initial_pairs_vec.size(); ++i){
        std::stack<IllumPair> node_stack;
        node_stack.push(initial_pairs_vec[i]);
        while(!node_stack.empty()){
            IllumPair curr = node_stack.top();
            node_stack.pop();

            if(!isIllumAware(curr, min_dist, error_thresh)){
                if(refineLTree(curr)){
                    if(curr.first->left != nullptr){
                        node_stack.push(std::make_pair(curr.first->left.get(), curr.second));
                    }

                    if(curr.first->right != nullptr){
                        node_stack.push(std::make_pair(curr.first->right.get(), curr.second));
                    }
                }
                else{
                    for(std::uint8_t i = 0; i < curr.second->children.size(); ++i){
                        if(curr.second->children[i] != nullptr){
                            node_stack.push(std::make_pair(curr.first, curr.second->children[i].get()));
                        }
                    }
                }
            }
            else{
                std::lock_guard<std::mutex> lock(illum_aware_mutex);
                illum_aware_pairs.push_back(curr);
            }
        }
    }

    return illum_aware_pairs;
}

std::uint64_t adaptiveVisibilitySampling(Scene* scene, LightTreeNode* light, OctreeNode<IllumcutSample>* curr_node, 
        std::unordered_map<std::uint32_t, bool>& visibility, std::mt19937& rng, float min_dist){
    if(curr_node->sample_indices.size() == 1){
        visibility[curr_node->sample_indices[0]] = sampleVisibility(scene, curr_node->representative().its, light->vpl, min_dist);
        return 1;
    }

    std::uint64_t taken_samples = 0;
    std::uint32_t max_samples = std::min(std::uint32_t(16), std::uint32_t(curr_node->sample_indices.size()));

    std::vector<float> dist_vals(8, 0.f);

    for(std::uint8_t i = 0; i < curr_node->children.size(); ++i){
        if(curr_node->children[i] != nullptr){
            dist_vals[i] = float(curr_node->children[i]->sample_indices.size());
        }
    }

    std::discrete_distribution<std::uint32_t> dist(dist_vals.begin(), dist_vals.end());

    std::vector<std::uint32_t> selected_indices;

    for(std::uint32_t i = 0; i < max_samples; ++i){
        std::uint32_t child_idx = dist(rng);
        std::uint32_t selected_idx = rand() % curr_node->children[child_idx]->sample_indices.size();
        std::uint32_t sample_idx = curr_node->children[child_idx]->sample_indices[selected_idx];
        selected_indices.push_back(sample_idx);
    }

    bool all_vis = true;
    bool not_all_invis = false;

    for(std::uint32_t i = 0; i < selected_indices.size(); ++i){
        bool vis;
        if(visibility.find(selected_indices[i]) != visibility.end()){
            vis = visibility[selected_indices[i]];
        }
        else{
            vis = sampleVisibility(scene, (*(curr_node->samples))[selected_indices[i]].its, light->vpl, min_dist);
            taken_samples++;
            visibility[selected_indices[i]] = vis;
        }

        all_vis = all_vis && vis;
        not_all_invis = not_all_invis || vis;
    }

    if(all_vis){
        for(std::uint32_t i = 0; i < curr_node->sample_indices.size(); ++i){
            visibility[curr_node->sample_indices[i]] = true;
        }
    }
    else if(!not_all_invis){
        for(std::uint32_t i = 0; i < curr_node->sample_indices.size(); ++i){
            visibility[curr_node->sample_indices[i]] = false;
        }
    }
    else{
        for(std::uint8_t i = 0; i < curr_node->children.size(); ++i){
            if(curr_node->children[i] != nullptr){
                taken_samples += adaptiveVisibilitySampling(scene, light, curr_node->children[i].get(), visibility, rng, min_dist);
            }
        }
    }

    return taken_samples;
}

//change to adaptive shadow sampling later
float renderIllumAwarePairs(const std::vector<IllumPair>& ilps, Scene* scene, float min_dist){
    Properties props("independent");
    ref<Sampler> sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
    sampler->configure();
    sampler->generate(Point2i(0));

    std::mutex render_mut;

    std::uint64_t total_samples = 0;
    std::uint64_t visibility_samples = 0;

    #pragma omp parallel for
    for(std::uint32_t i = 0; i < ilps.size(); ++i){
        std::unordered_map<std::uint32_t, bool> visibility;
        std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        std::uint64_t vis_samples = adaptiveVisibilitySampling(scene, ilps[i].first, ilps[i].second, visibility, rng, min_dist);

        {
            std::lock_guard<std::mutex> lock(render_mut);
            total_samples += ilps[i].second->sample_indices.size();
            visibility_samples += vis_samples;
        }

        for(std::uint32_t j = 0; j < ilps[i].second->sample_indices.size(); ++j){
            IllumcutSample& curr_sample = ilps[i].second->sample(j);
            std::uint32_t samples_taken;

            if(visibility[ilps[i].second->sample_indices[j]]){
                VPL vpl = ilps[i].first->vpl;
                vpl.P *= ilps[i].first->emission_scale;

                Spectrum col = sample(scene, sampler, curr_sample.its, curr_sample.ray, vpl, min_dist, 
                    false, 10, false, curr_sample.intersected_scene, true, false, samples_taken);

                {
                    std::lock_guard<std::mutex> lock(render_mut);
                    curr_sample.color += col;
                }
            }
        }
    }

    return float(visibility_samples) / total_samples;
}

void copySamplesToBuffer(std::uint8_t* output_image, const std::vector<IllumcutSample>& samples, Vector2i image_size,
    std::uint32_t spp){
    std::unordered_map<std::uint32_t, Spectrum> output;
    //#pragma omp parallel for
    for(std::uint32_t i = 0; i < samples.size(); ++i){
        std::uint32_t buffer_pos = samples[i].image_x + samples[i].image_y * image_size.x;
        if(output.find(buffer_pos) != output.end()){
            output[buffer_pos] += samples[i].color;
        }
        else output[buffer_pos] = samples[i].color;
    }

    for(auto iter = output.begin(); iter != output.end(); ++iter){
        Spectrum col = iter->second;
        col /= spp;

        float r, g, b;
        col.toSRGB(r, g, b);

        r = std::max(0.f, std::min(1.f, r));
        g = std::max(0.f, std::min(1.f, g));
        b = std::max(0.f, std::min(1.f, b));

        std::uint32_t buffer_pos = iter->first;

        output_image[buffer_pos*3] = r * 255 + 0.5f;
        output_image[buffer_pos*3 + 1] = g * 255 + 0.5f;
        output_image[buffer_pos*3 + 2] = b * 255 + 0.5f;
    }
}

bool IlluminationCutRenderer::render(Scene* scene, std::uint32_t spp, const RenderJob *job){
    srand(time(0));

    ref<Sensor> sensor = scene->getSensor();
    ref<Film> film = sensor->getFilm();

    auto size = film->getSize();
    if(size.x == 0 || size.y == 0){
        return true;
    }

    ref<Bitmap> output_bitmap = new Bitmap(Bitmap::ERGB, Bitmap::EUInt8, size);
    std::uint8_t* output_image = output_bitmap->getUInt8Data();
    memset(output_image, 0, output_bitmap->getBytesPerPixel() * size.x * size.y);

    std::unique_ptr<LightTree> light_tree(new LightTree(vpls_, min_dist_, 0, 0.f, false));
    std::cout << "Created light tree" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto receiver_root = constructOctree(scene, samples_, min_dist_, spp);
    std::cout << "Constructed octree" << std::endl;

    //computeUpperBounds(light_tree.get(), receiver_root.get(), scene, min_dist_, upper_distance_thresh_);
    computeUpperBounds2(light_tree.get(), receiver_root.get(), min_dist_, num_clusters_);
    std::cout << "Computed upper bounds" << std::endl;

    std::vector<IllumPair> illum_aware_pairs = getIlluminationAwarePairs(light_tree.get(), receiver_root.get(), min_dist_, error_threshold_);
    std::cout << "acquired " << illum_aware_pairs.size() << " illumination aware pairs" << std::endl;

    std::cout << "rendering..." << std::endl;
    float samplerate = renderIllumAwarePairs(illum_aware_pairs, scene, min_dist_);

    auto end = std::chrono::high_resolution_clock::now();

    std::vector<float> timings;
    std::vector<float> samplerates;
    timings.push_back(std::chrono::duration_cast<std::chrono::duration<float>>(end - start).count());
    samplerates.push_back(samplerate);

    writeOutputData(scene, "timings", false, timings, ',');
    writeOutputData(scene, "samplerates", false, samplerates, ',');

    std::cout << "Samplerate: " << samplerate << std::endl;

    /*std::stack<OctreeNode<IllumcutSample>*> node_stack;
    node_stack.push(receiver_root.get());
    while(!node_stack.empty()){
        OctreeNode<IllumcutSample>* curr = node_stack.top();
        node_stack.pop();

        if(curr->sample_indices.size() == 1){
            curr->representative().color = Spectrum(curr->upper_bound);
        }
        else{
            for(std::uint8_t i = 0; i < curr->children.size(); ++i){
                if(curr->children[i] != nullptr){
                    node_stack.push(curr->children[i].get());
                }
            }
        }
    }*/

    copySamplesToBuffer(output_image, samples_, size, spp);
    film->setBitmap(output_bitmap);

    return true;
}

MTS_NAMESPACE_END