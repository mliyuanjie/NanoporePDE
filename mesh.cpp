#define _USE_MATH_DEFINES

#include "mesh.h" 
#include "tools.h"
#include <algorithm>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <map>
#include <iomanip>

AtomNode::AtomNode(int id) {
	this->id = id;
}
AtomNode::~AtomNode() {
	if(this->next)
		delete this->next;
}

TreeNode::TreeNode(int level, double x, double y, double z, double dx) {
	this->level = level;
	this->x = x;
	this->y = y;
	this->z = z;
	this->dx = dx;
}
TreeNode::~TreeNode() {
	for (int i = 0; i < 8; i++) {
		if (this->child[i] != nullptr) {
			delete this->child[i];
			this->child[i] = nullptr;
		}
	}
	if (this->atomID != nullptr) {
		delete this->atomID;
		this->atomID = nullptr;
		this->atomEND = nullptr;
	}
	if (this->norm) {
		delete[] this->norm;
	}
	
}
OctaTree::OctaTree(int extradepth, int maxdepth, int mindepth, double x, double y, double z, double dx) {
	this->maxDepth = maxdepth;
	this->minDepth = mindepth;
	this->extraDepth = extradepth;
	this->minCell = dx;
	this->root = new TreeNode(0, x, y, z, dx / 2);
	this->root->type = NodeType::WaterVoids;
	this->initial(root);
}
OctaTree::~OctaTree() {
	delete this->root;
}
void OctaTree::initial(TreeNode* node) {
	int level = node->level;
	if (level >= minDepth) {
		//node->type = NodeType::Water;
		return;
	}
	this->refine(node);
	for (int i = 0; i < 8; i++) {
		this->initial(node->child[i]);
	}
};

void OctaTree::refine(TreeNode* node) {
	if (!node->isleaf())
		return;
	int level = node->level;
	double dx = node->dx;
	minCell = (dx < minCell) ? dx : minCell;
	this->size += 8;
	double x = node->x;
	double y = node->y;
	double z = node->z;
	for (int i = 0; i < 6; i++) {
		if (node->neighbor[i] && (node->level - node->neighbor[i]->level) > 0) {
			this->refine(node->neighbor[i]);
		}
	}
	node->child[0] = new TreeNode(level + 1, x - dx / 2, y - dx / 2, z - dx / 2, dx / 2);
	node->child[1] = new TreeNode(level + 1, x + dx / 2, y - dx / 2, z - dx / 2, dx / 2);
	node->child[2] = new TreeNode(level + 1, x - dx / 2, y + dx / 2, z - dx / 2, dx / 2);
	node->child[3] = new TreeNode(level + 1, x + dx / 2, y + dx / 2, z - dx / 2, dx / 2);
	node->child[4] = new TreeNode(level + 1, x - dx / 2, y - dx / 2, z + dx / 2, dx / 2);
	node->child[5] = new TreeNode(level + 1, x + dx / 2, y - dx / 2, z + dx / 2, dx / 2);
	node->child[6] = new TreeNode(level + 1, x - dx / 2, y + dx / 2, z + dx / 2, dx / 2);
	node->child[7] = new TreeNode(level + 1, x + dx / 2, y + dx / 2, z + dx / 2, dx / 2);
	for (int i = 0; i < 8; i++) {
		node->child[i]->parent = node;
		node->child[i]->type = node->type;
	}
	node->type = NodeType::External;
	//1
	node->child[0]->neighbor[0] = node->neighbor[0];
	if (node->neighbor[0] && !node->neighbor[0]->isleaf()) {
		node->child[0]->neighbor[0] = node->neighbor[0]->child[1];
		node->neighbor[0]->child[1]->neighbor[1] = node->child[0];
	}
	node->child[0]->neighbor[1] = node->child[1];
	node->child[0]->neighbor[2] = node->neighbor[2];
	if (node->neighbor[2] && !node->neighbor[2]->isleaf()) {
		node->child[0]->neighbor[2] = node->neighbor[2]->child[2];
		node->neighbor[2]->child[2]->neighbor[3] = node->child[0];
	}	
	node->child[0]->neighbor[3] = node->child[2];
	node->child[0]->neighbor[4] = node->neighbor[4];
	if (node->neighbor[4] && !node->neighbor[4]->isleaf()) {
		node->child[0]->neighbor[4] = node->neighbor[4]->child[4];
		node->neighbor[4]->child[4]->neighbor[5] = node->child[0];;
	}
	node->child[0]->neighbor[5] = node->child[4];
	//2
	node->child[1]->neighbor[0] = node->child[0];
	node->child[1]->neighbor[1] = node->neighbor[1];
	if (node->neighbor[1] && !node->neighbor[1]->isleaf()) {
		node->child[1]->neighbor[1] = node->neighbor[1]->child[0];
		node->neighbor[1]->child[0]->neighbor[0] = node->child[1];
	}
	node->child[1]->neighbor[2] = node->neighbor[2];
	if (node->neighbor[2] && !node->neighbor[2]->isleaf()) {
		node->child[1]->neighbor[2] = node->neighbor[2]->child[3];
		node->neighbor[2]->child[3]->neighbor[3] = node->child[1];
	}
	node->child[1]->neighbor[3] = node->child[3];
	node->child[1]->neighbor[4] = node->neighbor[4];
	if (node->neighbor[4] && !node->neighbor[4]->isleaf()) {
		node->child[1]->neighbor[4] = node->neighbor[4]->child[5];
		node->neighbor[4]->child[5]->neighbor[5] = node->child[1];
	}
	node->child[1]->neighbor[5] = node->child[5];
	//3
	node->child[2]->neighbor[0] = node->neighbor[0];
	if (node->neighbor[0] && !node->neighbor[0]->isleaf()) {
		node->child[2]->neighbor[0] = node->neighbor[0]->child[3];
		node->neighbor[0]->child[3]->neighbor[1] = node->child[2];
	}
	node->child[2]->neighbor[1] = node->child[3];
	node->child[2]->neighbor[2] = node->child[0];
	node->child[2]->neighbor[3] = node->neighbor[3];
	if (node->neighbor[3] && !node->neighbor[3]->isleaf()) {
		node->child[2]->neighbor[3] = node->neighbor[3]->child[0];
		node->neighbor[3]->child[0]->neighbor[2] = node->child[2];
	}
	node->child[2]->neighbor[4] = node->neighbor[4];
	if (node->neighbor[4] && !node->neighbor[4]->isleaf()) {
		node->child[2]->neighbor[4] = node->neighbor[4]->child[6];
		node->neighbor[4]->child[6]->neighbor[5] = node->child[2];
	}
	node->child[2]->neighbor[5] = node->child[6];
	//4
	node->child[3]->neighbor[0] = node->child[2];
	node->child[3]->neighbor[1] = node->neighbor[1];
	if (node->neighbor[1] && !node->neighbor[1]->isleaf()) {
		node->child[3]->neighbor[1] = node->neighbor[1]->child[2];
		node->neighbor[1]->child[2]->neighbor[0] = node->child[3];
	}
	node->child[3]->neighbor[2] = node->child[1];
	node->child[3]->neighbor[3] = node->neighbor[3];
	if (node->neighbor[3] && !node->neighbor[3]->isleaf()) {
		node->child[3]->neighbor[3] = node->neighbor[3]->child[1];
		node->neighbor[3]->child[1]->neighbor[2] = node->child[3];
	}
	node->child[3]->neighbor[4] = node->neighbor[4];
	if (node->neighbor[4] && !node->neighbor[4]->isleaf()) {
		node->child[3]->neighbor[4] = node->neighbor[4]->child[7];
		node->neighbor[4]->child[7]->neighbor[5] = node->child[3];
	}
	node->child[3]->neighbor[5] = node->child[7];
	//5
	node->child[4]->neighbor[0] = node->neighbor[0];
	if (node->neighbor[0] && !node->neighbor[0]->isleaf()) {
		node->child[4]->neighbor[0] = node->neighbor[0]->child[5];
		node->neighbor[0]->child[5]->neighbor[1] = node->child[4];
	}
	node->child[4]->neighbor[1] = node->child[5];
	node->child[4]->neighbor[2] = node->neighbor[2];
	if (node->neighbor[2] && !node->neighbor[2]->isleaf()) {
		node->child[4]->neighbor[2] = node->neighbor[2]->child[6];
		node->neighbor[2]->child[6]->neighbor[3] = node->child[4];
	}
	node->child[4]->neighbor[3] = node->child[6];
	node->child[4]->neighbor[4] = node->child[0];
	node->child[4]->neighbor[5] = node->neighbor[5];
	if (node->neighbor[5] && !node->neighbor[5]->isleaf()) {
		node->child[4]->neighbor[5] = node->neighbor[5]->child[0];
		node->neighbor[5]->child[0]->neighbor[4] = node->child[4];
	}
	//6
	node->child[5]->neighbor[0] = node->child[4];
	node->child[5]->neighbor[1] = node->neighbor[1];
	if (node->neighbor[1] && !node->neighbor[1]->isleaf()) {
		node->child[5]->neighbor[1] = node->neighbor[1]->child[4];
		node->neighbor[1]->child[4]->neighbor[0] = node->child[5];
	}
	node->child[5]->neighbor[2] = node->neighbor[2];
	if (node->neighbor[2] && !node->neighbor[2]->isleaf()) {
		node->child[5]->neighbor[2] = node->neighbor[2]->child[7];
		node->neighbor[2]->child[7]->neighbor[3] = node->child[5];
	}
	node->child[5]->neighbor[3] = node->child[7];
	node->child[5]->neighbor[4] = node->child[1];
	node->child[5]->neighbor[5] = node->neighbor[5];
	if (node->neighbor[5] && !node->neighbor[5]->isleaf()) {
		node->child[5]->neighbor[5] = node->neighbor[5]->child[1];
		node->neighbor[5]->child[1]->neighbor[4] = node->child[5];
	}
	//7
	node->child[6]->neighbor[0] = node->neighbor[0];
	if (node->neighbor[0] && !node->neighbor[0]->isleaf()) {
		node->child[6]->neighbor[0] = node->neighbor[0]->child[7];
		node->neighbor[0]->child[7]->neighbor[1] = node->child[6];
	}
	node->child[6]->neighbor[1] = node->child[7];
	node->child[6]->neighbor[2] = node->child[4];
	node->child[6]->neighbor[3] = node->neighbor[3];
	if (node->neighbor[3] && !node->neighbor[3]->isleaf()) {
		node->child[6]->neighbor[3] = node->neighbor[3]->child[4];
		node->neighbor[3]->child[4]->neighbor[2] = node->child[6];
	}
	node->child[6]->neighbor[4] = node->child[2];
	node->child[6]->neighbor[5] = node->neighbor[5];
	if (node->neighbor[5] && !node->neighbor[5]->isleaf()) {
		node->child[6]->neighbor[5] = node->neighbor[5]->child[2];
		node->neighbor[5]->child[2]->neighbor[4] = node->child[6];
	}
	//8
	node->child[7]->neighbor[0] = node->child[6];
	node->child[7]->neighbor[1] = node->neighbor[1];
	if (node->neighbor[1] && !node->neighbor[1]->isleaf()) {
		node->child[7]->neighbor[1] = node->neighbor[1]->child[6];
		node->neighbor[1]->child[6]->neighbor[0] = node->child[7];
	}
	node->child[7]->neighbor[2] = node->child[5];
	node->child[7]->neighbor[3] = node->neighbor[3];
	if (node->neighbor[3] && !node->neighbor[3]->isleaf()) {
		node->child[7]->neighbor[3] = node->neighbor[3]->child[5];
		node->neighbor[3]->child[5]->neighbor[2] = node->child[7];
	}
	node->child[7]->neighbor[4] = node->child[3];
	node->child[7]->neighbor[5] = node->neighbor[5];
	if (node->neighbor[5] && !node->neighbor[5]->isleaf()) {
		node->child[7]->neighbor[5] = node->neighbor[5]->child[3];
		node->neighbor[5]->child[3]->neighbor[4] = node->child[7];
	}
	
};
void OctaTree::_insertBall(TreeNode* node, double x0, double y0, double z0, double r0, int id0) {
	//find the max cell include the whole ball;
	double x = node->x;
	double y = node->y;
	double z = node->z;
	double dx = node->dx;
	double px = std::max(0.0, std::abs(x0 - x) - dx);
	double py = std::max(0.0, std::abs(y0 - y) - dx);
	double pz = std::max(0.0, std::abs(z0 - z) - dx);
	double d = sqrt(px * px + py * py + pz * pz);
	if (d > r0) {
		return;
	}
	//maintain the atom list in each node
	if (node->level >= extraDepth) {
		return;
	}
	if (node->isleaf()) { 
		this->refine(node); 
	}
	for (int i = 0; i < 8; i++) {
		this->_insertBall(node->child[i], x0, y0, z0, r0, id0);
	}
}
void OctaTree::_insertPore(TreeNode* node, double x0, double y0, double z0, double r0, double l0, double debye) {
	//find the max cell include the nanopore edge;
	//1. surface with large plane
	double x = node->x;
	double y = node->y;
	double z = node->z;
	double dx = node->dx;
	int cpos = cubicwithpore(x - x0, y - y0, z - z0, dx, r0, l0, debye);
	if ( cpos != 0) {
		return;
	}
	//maintain the atom list in each node
	if (node->level >= maxDepth) {
		return;
	}
	if (node->isleaf()) {
		this->refine(node);
	}
	for (int i = 0; i < 8; i++) {
		this->_insertPore(node->child[i], x0, y0, z0, r0, l0, debye);
	}
};

void OctaTree::visitSave(std::ofstream& f, std::ofstream& f2, TreeNode* node) {
	if (!node)
		return;
	if (node->isleaf()) {
		double x = node->x;
		double y = node->y;
		double z = node->z;
		double dx = node->dx;
		double pos[24] = {
			x - dx, y - dx, z - dx,
			x + dx, y - dx, z - dx,
			x - dx, y + dx, z - dx,
			x + dx, y + dx, z - dx,
			x - dx, y - dx, z + dx,
			x + dx, y - dx, z + dx,
			x - dx, y + dx, z + dx,
			x + dx, y + dx, z + dx,
		};
		f.write(reinterpret_cast<const char*>(pos), 24 * sizeof(double));
		double norm[3];
		if (node->norm) {
			norm[0] = node->norm[0];
			norm[1] = node->norm[1];
			norm[2] = node->norm[2];
		}
		else {
			norm[0] = 0;
			norm[1] = 0;
			norm[2] = 0;
		}
		double points[8] = {
			x, y, z, node->chargeDensity, double(node->type),
			norm[0], norm[1], norm[2]
		};
		f2.write(reinterpret_cast<const char*>(points), 8 * sizeof(double));
		return;
	}
	for (int i = 0; i < 8; i++) {
		this->visitSave(f, f2, node->child[i]);
	}
};

void OctaTree::save(std::string& fn, std::string& fn2) {
	std::ofstream outfile(fn, std::ios::binary);
	std::ofstream outfile2(fn2, std::ios::binary);
	//this->visitSave(outfile, outfile2, this->root);
	for (long long int i = 0; i < this->index.size(); i++) {
		TreeNode* node = this->index[i];
		double x = node->x;
		double y = node->y;
		double z = node->z;
		double dx = node->dx;
		double pos[24] = {
			x - dx, y - dx, z - dx,
			x + dx, y - dx, z - dx,
			x - dx, y + dx, z - dx,
			x + dx, y + dx, z - dx,
			x - dx, y - dx, z + dx,
			x + dx, y - dx, z + dx,
			x - dx, y + dx, z + dx,
			x + dx, y + dx, z + dx,
		};
		outfile.write(reinterpret_cast<const char*>(pos), 24 * sizeof(double));
		double norm[3];
		if (node->norm) {
			norm[0] = node->norm[0];
			norm[1] = node->norm[1];
			norm[2] = node->norm[2];
		}
		else {
			norm[0] = 0;
			norm[1] = 0;
			norm[2] = 0;
		}
		double points[8] = {
			x, y, z, node->chargeDensity, double(node->type),
			norm[0], norm[1], norm[2]
		};
		outfile2.write(reinterpret_cast<const char*>(points), 8 * sizeof(double));
	}
	outfile.close();
	outfile2.close();
}

void OctaTree::info() {
	std::vector<TreeNode*> stack;
	long long int n_ex = 0;
	long long int n_in = 0;
	stack.push_back(this->root);
	double maxcell = minCell;
	double porecell = 0;
	double protiencell = 0;
	long long int n_err_pore = 0;
	long long int n_err_protein = 0;
	long long int n_err = 0;
	long long int n_eff = 0;
	long long int n_sur_pore = 0;
	long long int n_sur_protein = 0;
	long long int n_bulk = 0;
	while (!stack.empty()) {
		TreeNode* node = stack.back();
		stack.pop_back();
		if (!node)
			continue;
		if (!node->isleaf()) {
			n_ex++;
		}
		else {
			n_in++;
		}
		if (node->isleaf()) {
			for (int i = 0; i < 6; i++) {
				if (!node->neighbor[i])
					continue;
				int level = node->neighbor[i]->level;
				if (abs(node->level - level) > 1) {
					n_err++;
				}
			}
			maxcell = (maxcell < 2 * node->dx) ? 2 * node->dx : maxcell;
		}
		if (node->type == NodeType::ProteinBoundary) {
			n_sur_protein++;
			if(protiencell == 0)
				protiencell = node->dx * 2;
			else if (protiencell != node->dx * 2) {
				n_err_protein++;
			}
			for (int i = 0; i < 6; i++) {
				if (node->neighbor[i]->type == NodeType::External) {
					n_err_protein++;
				}
			}
		}
		if (node->type == NodeType::PoreBoundary) {
			if (porecell == 0)
				porecell = node->dx * 2;
			else if (porecell != node->dx * 2) {
				n_err_pore++;
			}
			for (int i = 0; i < 6; i++) {
				if (node->neighbor[i] && node->neighbor[i]->type == NodeType::External) {
					n_err_pore++;
				}
			}
		}
		if (node->type != NodeType::External && node->type != NodeType::Pore) {
			n_eff++;
		}
		if (node->type == NodeType::Water) {
			n_bulk++;
		}
		if (node->type == NodeType::PoreBoundary)
			n_sur_pore++;
		for (int i = 0; i < 8; i++) {
			if (node->child[i]) {
				stack.push_back(node->child[i]);
			}
		}
	}
	std::cout << "INFO:----------------------" << std::endl;
	std::cout << "Leaf Node: " << n_in << ", Branch Node: " << n_ex << std::endl;
	std::cout << "min Cell Size: " << minCell << ", max Cell Size: " << maxcell << std::endl;
	std::cout << "Protein Cell Size: " << protiencell << ", Pore Cell Size: " << porecell << std::endl;
	std::cout << "Usage Total Node: " << n_eff << ", Usage Bulk Node: " << n_bulk <<
		std::endl << "Usage Pore Surface Node: " << n_sur_pore << ", Usage Protein Surface Node : " << n_sur_protein << std::endl;
	std::cout << "Error:---------------------" << std::endl;
	if (n_err_pore != 0) {
		std::cout << "Error Pore Node: " << n_err_pore << std::endl;
	}
	else if (n_err_protein != 0) {
		std::cout << "Error Protein Node: " << n_err_protein << std::endl;
	}
	else if (n_err != 0) {
		std::cout << "Error Node Level: " << n_err << std::endl;
	}
	else {
		std::cout << "Node Quality Good" << std::endl;
	}
	std::cout << "---------------------------" << std::endl;
}

void OctaTree::SAS(std::string& fn, double probe, double r, double l, double debye) {
	PNPNS::Atoms atoms;
	readPQR(atoms, fn, probe, 0, 0, 0);
	this->addProtein(atoms);
	this->addPore(r, l, debye);
	long long int n_sur_protien = 0;
	long long int n_sur_pore = 0;
	std::vector<TreeNode*> stack;
	TreeNode* node = this->root->child[6];
	while (!node->isleaf()) {
		node = node->child[5];
	}
	stack.push_back(node);
	while (!stack.empty()) {
		node = stack.back();
		stack.pop_back();
		if (!node)
			continue;
		if (node->type == NodeType::WaterVoids) {
			node->type = NodeType::Water;
			for (int i = 0; i < 6; i++) {
				if (!node->neighbor[i])
					continue;
				if (node->neighbor[i]->type < 0) {
					stack.push_back(node->neighbor[i]);
				}
				else if (!node->neighbor[i]->isleaf()) {
					for (int j = 0; j < 8; j++) {
						if (!node->neighbor[i]->child[j]->isleaf())
							continue;
						for (int k = 0; k < 6; k++) {
							if (node->neighbor[i]->child[j]->neighbor[k] == node &&
								node->neighbor[i]->child[j]->type < 0) {
								stack.push_back(node->neighbor[i]->child[j]);
							};
						}
							
					}
				}
			}
		}
		else if (node->type == NodeType::Pore) {
			node->type = NodeType::PoreBoundary;
			n_sur_pore++;
			for (int i = 0; i < 6; i++) {
				if (!node->neighbor[i])
					continue;
				if (node->neighbor[i]->type == NodeType::WaterVoids)
					stack.push_back(node->neighbor[i]);
				else if (!node->neighbor[i]->isleaf()) {
					for (int j = 0; j < 8; j++) {
						if (!node->neighbor[i]->child[j]->isleaf())
							continue;
						for (int k = 0; k < 6; k++) {
							if (node->neighbor[i]->child[j]->neighbor[k] == node &&
								node->neighbor[i]->child[j]->type == NodeType::WaterVoids) {
								stack.push_back(node->neighbor[i]->child[j]);
							};
						}

					}
				}
			}
		}
		else if (node->type == NodeType::Protein) {
			node->type = NodeType::ProteinBoundary;
			n_sur_protien++;
			for (int i = 0; i < 6; i++) {
				if (!node->neighbor[i])
					continue;
				if (node->neighbor[i]->type == NodeType::WaterVoids)
					stack.push_back(node->neighbor[i]);
				else if (!node->neighbor[i]->isleaf()) {
					for (int j = 0; j < 8; j++) {
						if (!node->neighbor[i]->child[j]->isleaf())
							continue;
						for (int k = 0; k < 6; k++) {
							if (node->neighbor[i]->child[j]->neighbor[k] == node &&
								node->neighbor[i]->child[j]->type == NodeType::WaterVoids) {
								stack.push_back(node->neighbor[i]->child[j]);
							};
						}

					}
				}
			}
		}
	}
	std::cout << "Protein Surface Node: " << n_sur_protien << ", Pore Surface Node: " << n_sur_pore << std::endl;
	this->_insertBoundary();
	this->_adjustBoundary();
	this->_setnorm(atoms, 0, 0, 0, r, l);
}

void OctaTree::addPore(double r, double l, double debye) {
	this->_insertPore(this->root, 0, 0, 0, r, l, debye);

	std::vector<TreeNode*> stack;
	TreeNode* node = this->root;
	stack.push_back(node);
	while (!stack.empty()) {
		node = stack.back();
		stack.pop_back();
		if (node->isleaf()) {
			double x = node->x;
			double y = node->y;
			double z = node->z;
			double dx = node->dx;
			if (pointnanopore(x, y, abs(z), r, l / 2) <= 0) {
				node->type = NodeType::Pore;
			}
		}
		else {
			for (int i = 0; i < 8; i++) {
				stack.push_back(node->child[i]);
			}
		}
	}
};
TreeNode* OctaTree::search(double x0, double y0, double z0) {
	TreeNode* node = this->root;
	double x;
	double y;
	double z;
	double dx;
	bool stop = false;
	while (!node->isleaf()) {
		for (int i = 0; i < 8; i++) {
			x = node->child[i]->x;
			y = node->child[i]->y;
			z = node->child[i]->z;
			dx = node->child[i]->dx;
			if (x0 >= (x - dx) && x0 <= (x + dx) &&
				y0 >= (y - dx) && y0 <= (y + dx) &&
				z0 >= (z - dx) && z0 <= (z + dx)) {
				node = node->child[i];
				break;
			}
		}
	}
	return node;
};

void OctaTree::addProtein(PNPNS::Atoms& atoms) {
	for (int i = 0; i < atoms.x.size(); i++) {
		double r = atoms.radius[i];
		double atomx = atoms.x[i];
		double atomy = atoms.y[i];
		double atomz = atoms.z[i];
		this->_insertBall(this->root, atomx, atomy, atomz, 3 * r, i);
	}
	std::vector<TreeNode*> stack;
	double invsqrt_pi = sqrt(0.5 / M_PI);
	for (int i = 0; i < atoms.x.size(); i++) {
		//std::cout << i << std::endl;
		double r = atoms.radius[i];
		double atomx = atoms.x[i];
		double atomy = atoms.y[i];
		double atomz = atoms.z[i];
		double charge = atoms.charge[i] * pow(invsqrt_pi / r, 3);
		TreeNode* node = this->search(atomx, atomy, atomz);
		stack.push_back(node);
		while (!stack.empty()) {
			node = stack.back();
			stack.pop_back();
			if (node && node->isleaf()) {
				double x = node->x;
				double y = node->y;
				double z = node->z;
				double d = sqrt((x - atomx) * (x - atomx) + (y - atomy) * (y - atomy) + (z - atomz) * (z - atomz));
				if (d <= r) {
					node->chargeDensity += charge * exp(-1.0 * d * d / r / r / 2.0);
					node->type = NodeType::ProteinBoundary;
					if (node->atomID == nullptr) {
						node->atomID = new AtomNode(i);
						node->atomEND = node->atomID;
					}
					else {
						node->atomEND->next = new AtomNode(i);
						node->atomEND = node->atomEND->next;
					}
					for (int i = 0; i < 6; i++) {
						if(node->neighbor[i]->type != NodeType::ProteinBoundary)
							stack.push_back(node->neighbor[i]);
					}
				}
				else if (d <= 3 * r) {
					node->chargeDensity += charge * exp(-1.0 * d * d / r / r / 2.0);
					node->type = NodeType::Water;
					for (int i = 0; i < 6; i++) {
						if (node->neighbor[i]->type != NodeType::ProteinBoundary && node->type != NodeType::Water)
							stack.push_back(node->neighbor[i]);
					}
				}
			}
		}
	}
	TreeNode* node = this->root;
	stack.push_back(node);
	while (!stack.empty()) {
		node = stack.back();
		stack.pop_back();
		if (!node)
			continue;
		if (node && node->isleaf()) {
			if(node->type == NodeType::Water)
				node->type = NodeType::WaterVoids;
			else if(node->type == NodeType::ProteinBoundary)
				node->type = NodeType::Protein;
		}
		else {
			for (int i = 0; i < 8; i++) {
				stack.push_back(node->child[i]);
			}
		}
	}
};

void OctaTree::_insertBoundary() {
	std::vector<TreeNode*> stack;
	TreeNode* node = this->root;
	stack.push_back(node);
	while (!stack.empty()) {
		node = stack.back();
		stack.pop_back();
		if (!node)
			continue;
		if (node->isleaf()) {
			if (!node->neighbor[4] || !node->neighbor[5]) {
				node->type = NodeType::BoundaryD;
			}
			else if (!(
				node->neighbor[0] && node->neighbor[1]
				&& node->neighbor[2] && node->neighbor[3]
				&& node->neighbor[4] && node->neighbor[5]
				) && node->type != NodeType::Pore) {
				node->type = NodeType::BoundaryN;
			}
		}
		else {
			for (int i = 0; i < 8; i++) {
				stack.push_back(node->child[i]);
			}
		}
	}
};

void OctaTree::_adjustBoundary() {
	long long int n_protein = 0;
	long long int n_pore = 0;
	std::vector<TreeNode*> stack;
	TreeNode* node = this->root;
	stack.push_back(node);
	while (!stack.empty()) {
		node = stack.back();
		stack.pop_back();
		if (!node)
			continue;
		if (node->isleaf() && (node->type == NodeType::ProteinBoundary)) {
			bool good = false;
			for (int i = 0; i < 6; i++) {
				if (node->neighbor[i]->type == NodeType::External) {
					std::cout << "Error Protein Surface Neignbour" << std::endl;
					return;
				}
				if (node->neighbor[i]->type == NodeType::Protein) {
					good = true;
				}
			}
			if (!good) {
				n_protein++;
				node->type = NodeType::Water;
			}
		}
		else if (node->isleaf() && (node->type == NodeType::PoreBoundary)) {
			bool good = false;
			for (int i = 0; i < 6; i++) {
				if (node->neighbor[i]->type == NodeType::External) {
					std::cout << "Error Pore Surface Neignbour" << std::endl;
					return;
				}
				if (node->neighbor[i]->type == NodeType::Pore) {
					good = true;
				}
			}
			if (!good) {
				n_pore++;
				node->type = NodeType::Water;
			}
		}
		else {
			for (int i = 0; i < 8; i++) {
				stack.push_back(node->child[i]);
			}
		}
	}
	std::cout << "Delete Protein Node: " << n_protein << ", Delete Pore Node: " << n_pore << std::endl;
}

void OctaTree::generateIndex() {
	std::queue<TreeNode*> stack;
	long long int id = 0;
	TreeNode* node = this->root->child[0];
	while (node->child[0]) {
		node = node->child[0];
	}
	stack.push(node);
	node->index = -2;
	while (!stack.empty()) {
		node = stack.front();
		stack.pop();
		if (!node)
			continue;
		if (node->isleaf()) {
			if (node->type != NodeType::Pore) {
				index.push_back(node);
				node->index = id;
				id++;
				for (int i = 0; i < 6; i++) {
					if (!node->neighbor[i])
						continue;
					if (node->neighbor[i]->isleaf() && node->neighbor[i]->index == -1) {
						node->neighbor[i]->index = -2;
						stack.push(node->neighbor[i]);
					}
					else if (!node->neighbor[i]->isleaf()) {
						for (int j = 0; j < 8; j++) {
							if (!node->neighbor[i]->child[j]->isleaf())
								continue;
							for (int k = 0; k < 6; k++) {
								if (node->neighbor[i]->child[j]->neighbor[k] == node &&
									node->neighbor[i]->child[j]->index == -1) {
									node->neighbor[i]->child[j]->index = -2;
									stack.push(node->neighbor[i]->child[j]);
								};
							}

						}
					}
				}
			}
		}
	}
	//std::cout << "number: " << id << std::endl;
}

void OctaTree::_setnorm(PNPNS::Atoms& atoms, double x0, double y0, double z0, double r, double l) {
	std::vector<TreeNode*> stack;
	TreeNode* node = this->root;
	stack.push_back(node);
	double rlarge2 = 1.2 * 1.2 * r * r;
	while (!stack.empty()) {
		node = stack.back();
		stack.pop_back();
		if (!node)
			continue;
		if (node->isleaf() && (node->type == NodeType::ProteinBoundary)) {
			double x = node->x;
			double y = node->y;
			double z = node->z;
			double* norm = new double[3];
			norm[0] = 0;
			norm[1] = 0;
			norm[2] = 0;
			AtomNode* atomi = node->atomID;
			while (atomi) {
				//double pre_coeff = atoms.charge[j] * ekbt * ee * 1e10 / (4.0 * M_PI * eps0);
				double r = atoms.radius[atomi->id];
				double atomx = atoms.x[atomi->id];
				double atomy = atoms.y[atomi->id];
				double atomz = atoms.z[atomi->id];
				double distan = sqrt((x - atomx) * (x - atomx) + (y - atomy) * (y - atomy) + (z - atomz) * (z - atomz));
				//double distan3 = pow(distan, 3);
				if (distan <= r) {
					norm[0] += (x - atomx);
					norm[1] += (y - atomy);
					norm[2] += (z - atomz);
				}
				atomi = atomi->next;
			}
			double r = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
			if (r > 0) {
				//double sign_t = (charge >= 0.0) ? 1 : -1;
				norm[0] /= r;
				norm[1] /= r;
				norm[2] /= r;
			}
			node->norm = norm;
		}
		else if (node->isleaf() && (node->type == NodeType::PoreBoundary)) {
			double x = node->x - x0;
			double y = node->y - y0;
			double z = node->z - z0;
			double d2 = x * x + y * y;
			double* norm = new double[3];
			norm[0] = 0;
			norm[1] = 0;
			norm[2] = 0;
			if (d2 >= rlarge2) {
				norm[2] = (z > 0) ? 1 : -1;
			}
			else if (d2 >= r * r && z <= l / 2 - 0.2 * r && z >= -l / 2 + 0.2 * r) {
				double distan = sqrt(d2);
				norm[0] = -x / distan;
				norm[1] = -y / distan;
			}
			else {
				if (z > 0) {
					 donutnormal(0.2 * r, l / 2 - 0.2 * r, 1.2 * r, x, y, z, norm);
				}
				else {
					donutnormal(0.2 * r, 0.2 * r - l / 2, 1.2 * r, x, y, z, norm);
				}
			}
			node->norm = norm;
		}
		else if (node->isleaf() && (node->type == NodeType::BoundaryN)) {
			double normtmp = 0;
			double* norm = new double[3];
			norm[0] = 0;
			norm[1] = 0;
			norm[2] = 0;
			if (node->neighbor[0]) {
				norm[0] += 1;
				normtmp += 1;
			}
			if (node->neighbor[1]) {
				norm[0] += -1;
				normtmp += 1;
			}
			if (node->neighbor[2]) {
				norm[1] += 1;
				normtmp += 1;
			}
			if (node->neighbor[3]) {
				norm[1] += -1;
				normtmp += 1;
			}
			if (node->neighbor[4]) {
				norm[2] += 1;
				normtmp += 1;
			}
			if (node->neighbor[5]) {
				norm[2] += -1;
				normtmp += 1;
			}
			normtmp = sqrt(normtmp);
			norm[0] /= normtmp;
			norm[1] /= normtmp;
			norm[2] /= normtmp;
			node->norm = norm;
		}
		else {
			for (int i = 0; i < 8; i++) {
				stack.push_back(node->child[i]);
			}
		}
	}
}

void OctaTree::check() {
	int graph[6][4] = {
		{ 1, 3, 5, 7 },
		{ 0, 2, 4, 6 },
		{ 2, 3, 6, 7 },
		{ 0, 1, 4, 5 },
		{ 4, 5, 6, 7 },
		{ 0, 1, 2, 3 }

	};
	long long int n = this->index.size();
	std::vector<TreeNode*>& nodes = this->index;
	std::vector<PNPNS::ETrip> Jlist;
	for (long long int i = 0; i < n; i++) {
		TreeNode* node = nodes[i];
		NodeType type = node->type;
		double res[3] = { 0, 0, 0 };
		double coeff[6] = { 1 };
		if (node->isleaf()) {
			for (int j = 0; j < 6; j++) {
				if (node->neighbor[j] && node->neighbor[j]->type == NodeType::Pore) {
					node->neighbor[j] = nullptr;
				}
				else if (node->neighbor[j] && node->neighbor[j]->type == NodeType::External) {
					for (int k = 0; k < 4; k++) {
						if (node->neighbor[j]->child[graph[j][k]]->type == NodeType::Pore) {
							node->neighbor[j] = nullptr;
						}
					}
				}
			}
		}
	}
}

template <typename T>
void SwapEnd(T& var)
{
	char* varArray = reinterpret_cast<char*>(&var);
	for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
		std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

void OctaTree::toVTK(std::string& fn) {
	using Point = std::tuple<float, float, float>;
	std::ofstream outFile(fn, std::ios::binary);
	outFile << "# vtk DataFile Version 2.0\n";
	outFile << "Protein in Nanopore\n";
	outFile << "BINARY\n";
	outFile << "DATASET UNSTRUCTURED_GRID\n";
	//build the cell and points vector;
	long long int n = this->index.size();
	long long int vertices_id = 0;
	std::vector<TreeNode*>& nodes = this->index;
	std::vector<float> vertices;
	std::vector<std::vector<long long int>> cells(n, std::vector<long long int> (9, 0));
	std::map<Point, long long int> uniquePoints;

	for (long long int i = 0; i < n; i++) {
		TreeNode* node = nodes[i];
		cells[node->index][0] = 8;
		float x = node->x;
		float y = node->y;
		float z = node->z;
		float dx = node->dx;
		std::vector<Point> pos = {
			{x - dx, y - dx, z - dx},
			{x + dx, y - dx, z - dx},
			{x - dx, y + dx, z - dx},
			{x + dx, y + dx, z - dx},
			{x - dx, y - dx, z + dx},
			{x + dx, y - dx, z + dx},
			{x - dx, y + dx, z + dx},
			{x + dx, y + dx, z + dx},
		};
		for (int j = 0; j < 8; j++) {
			if (uniquePoints.find(pos[j]) == uniquePoints.end()) {
				uniquePoints[pos[j]] = vertices_id;
				vertices.push_back(std::get<0>(pos[j]));
				vertices.push_back(std::get<1>(pos[j]));
				vertices.push_back(std::get<2>(pos[j]));
				cells[node->index][j + 1] = vertices_id;
				vertices_id++;
			}
			else {
				long long int vertices_id_old = uniquePoints[pos[j]];
				cells[node->index][j + 1] = vertices_id_old;
			}
		}
	}
	uniquePoints.clear();

	// Write vertices
	long long int numPoints = vertices.size() / 3;
	outFile << "POINTS " << numPoints << " float\n";
	for (auto& vertex : vertices) {
		SwapEnd(vertex);
		outFile.write(reinterpret_cast<char*>(&vertex), sizeof(vertex));
	}

	// Write cells
	uint32_t totalCellData = 9 * cells.size();
	outFile << "\nCELLS " << cells.size() << " " << totalCellData << "\n";
	for (const auto& cell : cells) {
		for (int i = 0; i < 9; i++) { 
			long long int vertex = cell[i];
			uint32_t vertexBE = static_cast<uint32_t>(vertex);
			SwapEnd(vertexBE);
			outFile.write(reinterpret_cast<char*>(&vertexBE), sizeof(vertexBE));
		}
	}

	outFile << "\nCELL_TYPES " << cells.size() << "\n";
	uint32_t cellType = 12;
	SwapEnd(cellType);
	for (long long int i = 0; i < cells.size(); ++i) {
		outFile.write(reinterpret_cast<char*>(&cellType), sizeof(cellType));
	}
	outFile.close();
	std::cout << "save to vtk" << std::endl;
}

