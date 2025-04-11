#pragma once

#include <string>
#include <vector>
#include "tools.h"


enum NodeType {
	External = 0,
	BoundaryD = 4,
	BoundaryN = -4,
	Water = 1,
	WaterVoids = -1,
	Pore = -2,
	PoreBoundary = 2,
	Protein = -3,
	ProteinBoundary = 3,
};
class AtomNode {
public:
	AtomNode(int id);
	~AtomNode();
	int id = -1;
	AtomNode* next = nullptr;
};
class TreeNode {
public:
	TreeNode(int level, double x, double y, double z, double dx);
	~TreeNode();
	bool isleaf() { return (this->child[0] == nullptr) ? true : false; };
	long long int index = -1;
	int level;
	NodeType type = NodeType::External;
	double x;
	double y;
	double z;
	double dx;
	double chargeDensity = 0.0;
	AtomNode* atomID = nullptr;
	AtomNode* atomEND = nullptr;
	TreeNode* child[8] = { nullptr };
	TreeNode* neighbor[6] = { nullptr };
	TreeNode* parent = nullptr;
	double* norm = nullptr;
};

class OctaTree {
public:
	OctaTree(int extradepth, int maxdepth, int mindepth, double x, double y, double z, double dx);
	~OctaTree();
	void initial(TreeNode* node);

	TreeNode* root = nullptr;
	std::vector<TreeNode*> index;
	int maxDepth;
	int minDepth;
	int extraDepth;
	double minCell;
	long long int size = 0;

	void refine(TreeNode* node);
	void addPore(double r, double l, double debye);
	void addProtein(PNPNS::Atoms& atoms);
	void _insertBoundary();
	void _insertBall(TreeNode* node, double x0, double y0, double z0, double r0, int id0);
	void _insertPore(TreeNode* node, double x0, double y0, double z0, double r0, double l0, double debye);
	void _adjustBoundary();
	TreeNode* search(double x0, double y0, double z0);
	void save(std::string& fn, std::string& fn2);
	void visitSave(std::ofstream& f, std::ofstream& f2, TreeNode* node);
	void generateIndex();
	void info();
	void SAS(std::string& fn, double probe, double r, double l, double debye);
	void _setnorm(PNPNS::Atoms& atoms, double x0, double y0, double z0, double r, double l);
	void check();
	void toVTK(std::string& fn);
};
