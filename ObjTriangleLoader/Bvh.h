#pragma once

#include<memory>
#include "Triangle.h"

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


};