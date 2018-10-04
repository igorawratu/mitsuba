#include "lighttree.h"

#include <tuple>
#include <queue>
#include <functional>

MTS_NAMESPACE_BEGIN

LightTree::LightTree() : point_tree_root_(nullptr), directional_tree_root_(nullptr), oriented_tree_root_(nullptr){

}

void calculateSimilarityMatrices(const std::vector<VPL>& point_vpls, const std::vector<VPL>& directional_vpls, const std::vector<VPL>& oriented_vpls,
	std::vector<float>& point_similarity_matrix, std::vector<float>& direction_similarity_matrix, std::vector<float>& oriented_similarity_matrix) {
	
	point_similarity_matrix.clear();
	if (point_vpls.size() > 0) {
		point_similarity_matrix.insert(point_similarity_matrix.begin(), point_vpls.size() * point_vpls.size(), 0.f);
	}

	direction_similarity_matrix.clear();
	if (directional_vpls.size() > 0) {
		direction_similarity_matrix.insert(direction_similarity_matrix.begin(), directional_vpls.size() * directional_vpls.size(), 0.f);
	}

	oriented_similarity_matrix.clear();
	if (oriented_vpls.size() > 0) {
		oriented_similarity_matrix.insert(oriented_similarity_matrix.begin(), oriented_vpls.size() * oriented_vpls.size(), 0.f);
	}

	for (std::uint32_t i = 0; i < point_vpls.size(); ++i) {
		for (std::uint32_t j = 0; j < point_vpls.size(); ++j) {
			std::uint32_t index = j + i * point_vpls.size();
			point_similarity_matrix[index] = (point_vpls[i].its.p - point_vpls[j].its.p).length();
		}
	}

	for (std::uint32_t i = 0; i < directional_vpls.size(); ++i) {
		for (std::uint32_t j = 0; j < directional_vpls.size(); ++j) {
			std::uint32_t index = j + i * directional_vpls.size();
			direction_similarity_matrix[index] = dot(directional_vpls[i].its.shFrame.n, directional_vpls[j].its.shFrame.n);
		}
	}

	for (std::uint32_t i = 0; i < oriented_vpls.size(); ++i) {
		for (std::uint32_t j = 0; j < oriented_vpls.size(); ++j) {
			std::uint32_t index = j + i * oriented_vpls.size();
			oriented_similarity_matrix[index] = (oriented_vpls[i].its.p - oriented_vpls[j].its.p).length() + dot(oriented_vpls[i].its.shFrame.n, oriented_vpls[j].its.shFrame.n);
		}
	}
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
		
	}

	typedef std::tuple<float, std::uint32_t, std::uint32_t> SimEntry;

	auto comparator = [](SimEntry l, SimEntry r){
		return std::get<0>(l) < std::get<0>(r);
	};

	std::priority_queue<SimEntry, std::vector<SimEntry>, decltype(comparator)> similarity_matrix;
	
	while(true){
		if(nodes[current_level].size() == 1){
			break;
		}

		for(int i = 0; i < nodes[current_level].size(); ++i){
			for(int j = i + 1; j < nodes[current_level].size(); ++j){
				float d = 0.f;
				if(vpl_type != EDirectionalEmitterVPL){
					d += (nodes[current_level][i]->br - nodes[current_level][i]->tl).length();
				}
				
				if(vpl_type != EPointEmitterVPL){
					d += (1.f - dot(nodes[current_level][i]->cone_ray1, nodes[current_level][i]->cone_ray2)) / 0.5f;
				}

				similarity_matrix.push(std::make_tuple(d, i, j));
			}
		}

		std::uint32_t nodes_left = nodes[current_level].size();
		std::uint8_t next_level = (current_level + 1) % 2;

		while(nodes_left > 0){
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
					dp2 = dot(nodes[next_level].back()->cone_ray1, c2->cone_ray1);
					dp3 = dot(nodes[next_level].back()->cone_ray1, c2->cone_ray2);
					
				}
				
			}
		}
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

	calculateSimilarityMatrices(point_vpls_, directional_vpls_, oriented_vpls_, point_similarity_matrix_, directional_similarity_matrix_, oriented_similarity_matrix_);
}

LightTree::LightTree(const LightTree& other) {

}

LightTree::LightTree(const LightTree&& other) {

				continue;
		}
	}

	calculateSimilarityMatrices(point_vpls_, directional_vpls_, oriented_vpls_, point_similarity_matrix_, directional_similarity_matrix_, oriented_similarity_matrix_);
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