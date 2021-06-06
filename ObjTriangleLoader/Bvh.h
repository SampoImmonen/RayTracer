#pragma once

#include <memory>
#include <fstream>
#include <string>
#include <iostream>

#include "RTTriangle.h"
#include "Utils.h"

namespace Helmi {

struct BvhNode {
	BoundingBox bb;
	uint16_t start, end;
	BvhNode* left;
	BvhNode* right;
	int splitAxis;

	BvhNode(): bb(), start(0), end(0), left(nullptr), right(nullptr){}
	~BvhNode() {
		delete left;
		delete right;
	}
};

struct FlatBvhNode {
	BoundingBox bb;
	union {
		int start;
		int childOffset;
	};
	int16_t num_triangles;
	uint8_t axis;
	uint8_t pad[1];
};


inline std::ofstream& operator<<(std::ofstream& of, const FlatBvhNode& n) {
	of << n.bb.min;
	of << n.bb.max;
	of << n.start << " " << n.num_triangles << " " << static_cast<unsigned>(n.axis) << std::endl;

	return of;
}

inline std::ifstream& operator>>(std::ifstream& is, FlatBvhNode& n) {
	is >> n.bb.min >> n.bb.max >> n.start >> n.num_triangles >> n.axis;
	return is;
}


//overload glm::vec3 for serialization

class Bvh {

public:

	Bvh(){};
	~Bvh() {
		delete rootNode;
	}

	Bvh& operator=(Bvh&& other) {
		std::swap(rootNode, other.rootNode);
		std::swap(indices, other.indices);
		return *this;
	}

	BvhNode& root() { return *rootNode; }
	const BvhNode& root() const { return *rootNode; }
	uint16_t	getIndex(uint16_t index) const { return indices[index]; }

	BvhNode* rootNode;
	std::vector<uint16_t> indices;
	int n_nodes = 0;
	std::vector<FlatBvhNode> flatnodes;

	void initFlatBvh() {
		flatnodes.resize(n_nodes, FlatBvhNode());
	}

	int FlattenBvhTree(BvhNode* N, int* offset) {

		FlatBvhNode* v = &flatnodes[*offset];
		v->bb = N->bb;
		int myOffset = (*offset)++;
		if (N->left == nullptr && N->right == nullptr) {
			//handle leaf nodes
			v->start = N->start;
			v->num_triangles = N->end - N->start;
		}
		else {
			//handle interior nodes
			//we dont care about this at this point (maybe later)
			v->axis = N->splitAxis;
			v->num_triangles = 0;
			FlattenBvhTree(N->left, offset);
			v->childOffset = FlattenBvhTree(N->right, offset);
		}
		return myOffset;
	}

	void save(const std::string& filepath) {

		std::cout << "saving Bvh of scene into: " + filepath << "\n";
		std::ofstream savefile(filepath.c_str());

		//store number of indices and number of nodes
		savefile << indices.size() << "\n";
		savefile << n_nodes << "\n";
		//save indices
		for (auto& i : indices) {
			savefile << i << " ";
		}
		//save FlatBvhNodes
		for (auto& n : flatnodes) {
			savefile << n;
		}
		
	}

	void load(const std::string& filepath) {
		

		std::cout << "loading bvh from path: " + filepath << "\n";
		std::ifstream loadfile(filepath.c_str());
		
		if (loadfile) {
			// load number of indices and nodes
			size_t num_indices;
			loadfile >> num_indices;
			loadfile >> n_nodes;

			// load indices
			std::cout << "loading indices\n";
			indices.resize(num_indices);
			for (size_t i = 0; i < num_indices; ++i) {
				loadfile >> indices[i];
			}

			// load FlatBvhNodes
			std::cout << "loading nodes\n";
			flatnodes.resize(n_nodes);
			for (int i = 0; i < n_nodes; ++i) {
				FlatBvhNode n;
				loadfile >> n;
				flatnodes[i] = n;
			}
		}
		else {
			std::cout << "Failed to open file " + filepath << "\n";
		}
	}
};

}