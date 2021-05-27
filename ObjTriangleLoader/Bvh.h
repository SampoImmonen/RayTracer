#pragma once

#include<memory>
#include "Triangle.h"

struct BvhNode {
	BoundingBox bb;
	uint16_t start, end;
	std::unique_ptr<BvhNode> left;
	std::unique_ptr<BvhNode> right;

	BvhNode(): bb(), start(0), end(0), left(nullptr), right(nullptr){}

};

struct FlatBvhNode {
	BoundingBox bb;

};

class Bvh {

public:

	Bvh(){};

	Bvh& operator=(Bvh&& other) {
		std::swap(rootNode, other.rootNode);
		std::swap(indices, other.indices);
		return *this;
	}

	BvhNode& root() { return *rootNode; }
	const BvhNode& root() const { return *rootNode; }
	uint16_t	getIndex(uint16_t index) const { return indices[index]; }

	std::unique_ptr<BvhNode> rootNode;
	std::vector<uint16_t> indices;
	int n_nodes = 0;
};