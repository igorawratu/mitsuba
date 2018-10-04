#include "lighttree.h"

#include <tuple>
#include <queue>
#include <functional>

MTS_NAMESPACE_BEGIN

LightTree::LightTree() : point_tree_root_(nullptr), directional_tree_root_(nullptr), oriented_tree_root_(nullptr){

}

std::unique_ptr<LightTreeNode> createLightTree(const std::vector<VPL>& vpls, EVPLType vpl_type) {
	if(vpls.size() == 0 || vpls[0].type != vpl_type){
		return nullptr;
	}

	std::array<std::vector<std::unique_ptr<LightTreeNode>>, 2> nodes;
	std::uint8_t current_level = 0;
	
	for(int i = 0; i < vpls.size(); ++i){
		nodes[current_level].emplace_back(new LightTreeNode());

		if(vpl_type != EDirectionalEmitterVPL){
			nodes[current_level][i]->tl = vpls[i].its.p;
			nodes[current_level][i]->br = vpls[i].its.p;
		}
		
		if(vpl_type != EPointEmitterVPL){
			nodes[current_level][i]->cone_ray1 = vpls[i].its.shFrame.n;
			nodes[current_level][i]->cone_ray2 = vpls[i].its.shFrame.n;
		}
		
		nodes[current_level][i]->vpl = vpls[i];
		nodes[current_level][i]->emission_scale = 1.f;
	}

	typedef std::tuple<float, std::uint32_t, std::uint32_t> SimEntry;

	auto comparator = [](SimEntry l, SimEntry r){
		return std::get<0>(l) < std::get<0>(r);
	};

	Properties props("halton");
	props.setInteger("scramble", 0);
	Sampler *sampler = static_cast<Sampler*>(PluginManager::getInstance()->createObject(MTS_CLASS(Sampler), props));
	sampler->configure();
	sampler->generate(Point2i(0));

	while(true){
		if(nodes[current_level].size() == 1){
			break;
		}

		std::priority_queue<SimEntry, std::vector<SimEntry>, decltype(comparator)> similarity_matrix;

		for(int i = 0; i < nodes[current_level].size(); ++i){
			for(int j = i + 1; j < nodes[current_level].size(); ++j){
				float d = 0.f;
				if(vpl_type != EDirectionalEmitterVPL){
					d += (nodes[current_level][i]->br - nodes[current_level][i]->tl).length();
				}
				
				if(vpl_type != EPointEmitterVPL){
					d += (1.f - dot(nodes[current_level][i]->cone_ray1, nodes[current_level][i]->cone_ray2));
				}

				similarity_matrix.push(std::make_tuple(d, i, j));
			}
		}

		std::uint32_t nodes_left = nodes[current_level].size();
		std::uint8_t next_level = (current_level + 1) % 2;

		while(nodes_left > 0){
			if (nodes_left == 1) {
				for (std::uint32_t i = 0; i < nodes[current_level].size(); ++i) {
					if (nodes[current_level][i] != nullptr) {
						nodes[next_level].push_back(std::move(nodes[current_level][i]));
						nodes[current_level][i] = nullptr;
						nodes_left--;
						break;
					}
				}

				continue;
			}

			SimEntry entry = similarity_matrix.top();
			similarity_matrix.pop();

			std::uint32_t id1 = std::get<1>(entry), id2 = std::get<2>(entry);

			if(nodes[current_level][id1] != nullptr && 
				nodes[current_level][id2] != nullptr){
				nodes[next_level].emplace_back(new LightTreeNode());
				LightTreeNode *c1 = nodes[current_level][id1].get(), *c2 = nodes[current_level][id2].get();
				
				if(vpl_type != EDirectionalEmitterVPL){
					nodes[next_level].back()->tl = Point(std::max(c1->tl[0], c2->tl[0]), 
						std::max(c1->tl[1], c2->tl[1]), 
						std::max(c1->tl[2], c2->tl[2]));
					nodes[next_level].back()->br = Point(std::min(c1->br[0], c2->br[0]), 
						std::min(c1->br[1], c2->br[1]), 
						std::min(c1->br[2], c2->br[2]));
				}
				
				//need to ideally find a better method to union 2 cones
				if(vpl_type != EPointEmitterVPL){
					nodes[next_level].back()->cone_ray1 = c1->cone_ray1;
					nodes[next_level].back()->cone_ray2 = c1->cone_ray2;
					
					float dp1 = dot(nodes[next_level].back()->cone_ray1, nodes[next_level].back()->cone_ray2);
					float dp2 = dot(nodes[next_level].back()->cone_ray1, c2->cone_ray1);
					float dp3 = dot(nodes[next_level].back()->cone_ray1, c2->cone_ray2);

					if(dp2 > dp1 || dp3 > dp1){
						nodes[next_level].back()->cone_ray2 = d2 > dp3 ? c2->cone_ray1 : c2->cone_ray2;
					}

					dp1 = dot(nodes[next_level].back()->cone_ray1, nodes[next_level].back()->cone_ray2);
					dp2 = dot(nodes[next_level].back()->cone_ray2, c2->cone_ray1);
					dp3 = dot(nodes[next_level].back()->cone_ray2, c2->cone_ray2);

					if (dp2 > dp1 || dp3 > dp1) {
						nodes[next_level].back()->cone_ray1 = dp2 > dp3 ? c2->cone_ray1 : c2->cone_ray2;
					}
				}
				
				float c1_intensity = c1->vpl.P.getLuminance() * c1->emission_scale;
				float c2_intensity = c2->vpl.P.getLuminance() * c2->emission_scale;

				float sample = sampler->next1D();

				if (sample < c1_intensity / (c1_intensity + c2_intensity)) {
					nodes[next_level].back().vpl = c1->vpl;
					nodes[next_level].back().emission_scale = c1->emission_scale + c2_intensity / c1->vpl.P.getLuminance();
				}
				else {
					nodes[next_level].back().vpl = c2->vpl;
					nodes[next_level].back().emission_scale = c2->emission_scale + c1_intensity / c2->vpl.P.getLuminance();
				}

				nodes[next_level].back().left = std::move(nodes[current_level][id1]);
				nodes[next_level].back().right = std::move(nodes[current_level][id2]);

				nodes[current_level][id1] = nullptr;
				nodes[current_level][id2] = nullptr;

				nodes_left -= 2;
			}
		}

		current_level = next_level;
	}

	return std::move(nodes[current_level][0]);
}

LightTree::LightTree(const std::vector<VPL>& vpls) {
	for (std::uint32_t i = 0; i < vpls.size(); ++i) {
		switch (vpls[i].type) {
			case EPointEmitterVPL:
				point_vpls_.push_back(vpls[i]);
				break;
			case ESurfaceVPL:
				oriented_vpls_.push_back(vpls[i]);
				break;
			case EDirectionalEmitterVPL:
				directional_vpls_.push_back(vpls[i]);
				break;
			default:
				continue;
		}
	}

	
}

LightTree::LightTree(const LightTree& other) {

}

LightTree::LightTree(const LightTree&& other) {

}

LightTree::LightTree(const LightTree& other) {

}

LightTree::LightTree(const LightTree&& other) {

}

LightTree& LightTree::operator = (const LightTree& other) {

}

LightTree& LightTree::operator = (LightTree&& other) {

}
}

LightTree& LightTree::operator = (const LightTree& other) {

}

LightTree& LightTree::operator = (LightTree&& other) {

}

LightTree::~LightTree() {

}

void LightTree::setVPLs(const std::vector<VPL>& vpls) {

}

std::vector<VPL> LightTree::getClusteringForPoint(Intersection its) {

}

MTS_NAMESPACE_ENDvectorvector