#include "operator.h" 

PDEOperator::PDEOperator(OctaTree* mesh) {
	this->mesh = mesh;
	long long int n = mesh->index.size();
	//this->JMatrix = new PartialList[n];
}

void PDEOperator::setPhysics(double conc, double temperature, int z0, int z1, double Diff0, double Diff1, double rho, double viscosity) {
	this->kbt = 1.380649e-23 * (273.15 + temperature);
	this->z0 = z0;
	this->z1 = z1;
	this->diffu0 = Diff0;
	this->diffu1 = Diff1;
	this->rho = rho;
	this->viscosity = viscosity;
	this->c0 = conc;
	this->ekbt = ee / kbt;

	//possion const term for charge density
	double charge_coeff = ee * avogadro * ekbt / eps0 * 1e-20 * c0;
	this->charge_coeff0 = charge_coeff * z0;
	this->charge_coeff1 = charge_coeff * z1;
	//nernst planck couple fluid
	this->np_ns_coeff0 = 1e-14 / diffu0;
	this->np_ns_coeff1 = 1e-14 / diffu1;
	//navier stock equation for viscosity , 
	this->Re_inv = viscosity / 1e-6 / P0;
	// navier coupled with space force
	this->ns_pnp_coeff = eps0 / ekbt / P0 * 1e20 / ekbt;
};
void PDEOperator::grad(TreeNode* node, double* x, double* y) {
	double x1 = x[node->index];
	for (int i = 0; i < 3; i++) {
		double x0 = 0;
		double x2 = 0;
		if (node->neighbor[i * 2] && node->neighbor[i * 2 + 1]) {
			double dx0;
			double dx1;
			if (node->neighbor[i * 2]->isleaf()) {
				x0 = x[node->neighbor[i * 2]->index];
				dx0 = node->dx + node->neighbor[i * 2]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x0 += x[node->neighbor[i * 2]->child[graph[i * 2][j]]->index];
				}
				x0 /= 4;
				dx0 = 1.5 * node->dx;
			}
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				x2 = x[node->neighbor[i * 2 + 1]->index];
				dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x2 += x[node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index];
				}
				x2 /= 4;
				dx1 = 1.5 * node->dx;
			}
			y[i] = (dx0 * dx0 * x2 + (dx1 * dx1 - dx0 * dx0) * x1 - dx1 * dx1 * x0) / (dx0 * dx1 * (dx0 + dx1));
		}
		else if (node->neighbor[i * 2]) {
			double dx;
			if (node->neighbor[i * 2]->isleaf()) {
				x0 = x[node->neighbor[i * 2]->index];
				dx = node->dx + node->neighbor[i * 2]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x0 += x[node->neighbor[i * 2]->child[graph[i * 2][j]]->index];
				}
				x0 /= 4;
				dx = 1.5 * node->dx;
			}
			y[i] = (x1 - x0) / dx;
		}
		else if (node->neighbor[i * 2 + 1]) {
			double dx;
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				x2 = x[node->neighbor[i * 2 + 1]->index];
				dx = node->dx + node->neighbor[i * 2 + 1]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x2 += x[node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index];
				}
				x2 /= 4;
				dx = 1.5 * node->dx;
			}
			y[i] = (x2 - x1) / dx;
		}
	}
	return;
}
void PDEOperator::lapl(TreeNode* node, double* x, double* y, double* coeff) {
	double x1 = x[node->index];
	for (int i = 0; i < 3; i++) {
		double x0 = 0;
		double x2 = 0;
		if (node->neighbor[i * 2] && node->neighbor[i * 2 + 1]) {
			double dx0;
			double dx1;
			if (node->neighbor[i * 2]->isleaf()) {
				x0 = x[node->neighbor[i * 2]->index];
				dx0 = node->dx + node->neighbor[i * 2]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x0 += x[node->neighbor[i * 2]->child[graph[i * 2][j]]->index];
				}
				x0 /= 4;
				dx0 = 1.5 * node->dx;
			}
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				x2 = x[node->neighbor[i * 2 + 1]->index];
				dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x2 += x[node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index];
				}
				x2 /= 4;
				dx1 = 1.5 * node->dx;
			}
			y[i] = 2 * ((x2 - x1) / dx1 * coeff[i * 2 + 1] - (x1 - x0) / dx0 * coeff[i * 2]) / (dx0 + dx1);
		}
		else if (node->neighbor[i * 2]) {
			double dx;
			if (node->neighbor[i * 2]->isleaf()) {
				x0 = x[node->neighbor[i * 2]->index];
				dx = node->dx + node->neighbor[i * 2]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x0 += x[node->neighbor[i * 2]->child[graph[i * 2][j]]->index];
				}
				x0 /= 4;
				dx = 1.5 * node->dx;
			}
			y[i] = (x0 - x1) / dx / dx * coeff[i * 2];
		}
		else if (node->neighbor[i * 2 + 1]) {
			double dx;
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				x2 = x[node->neighbor[i * 2 + 1]->index];
				dx = node->dx + node->neighbor[i * 2 + 1]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x2 += x[node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index];
				}
				x2 /= 4;
				dx = 1.5 * node->dx;
			}
			y[i] = (x2 - x1) / dx / dx * coeff[i * 2 + 1];
		}
	}
	return;
}
void PDEOperator::grad_direction(TreeNode* node, double* x, double* y, double* direction) {
	double x1 = x[node->index];
	for (int i = 0; i < 3; i++) {
		double x0 = 0;
		double x2 = 0;
		if (direction[i] == 0) {
			double dx0 = 0;
			double dx1 = 0;
			if (node->neighbor[i * 2] && node->neighbor[i * 2]->isleaf()) {
				x0 = x[node->neighbor[i * 2]->index];
				dx0 = node->dx + node->neighbor[i * 2]->dx;
			}
			else if(node->neighbor[i * 2]) {
				for (int j = 0; j < 4; j++) {
					x0 += x[node->neighbor[i * 2]->child[graph[i * 2][j]]->index];
				}
				x0 /= 4;
				dx0 = 1.5 * node->dx;
			}
			if (node->neighbor[i * 2 + 1] && node->neighbor[i * 2 + 1]->isleaf()) {
				x2 = x[node->neighbor[i * 2 + 1]->index];
				dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
			}
			else if(node->neighbor[i * 2 + 1]) {
				for (int j = 0; j < 4; j++) {
					x2 += x[node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index];
				}
				x2 /= 4;
				dx1 = 1.5 * node->dx;
			}
			y[i] = (dx0 * dx0 * x2 + (dx1 * dx1 - dx0 * dx0) * x1 - dx1 * dx1 * x0) / (dx0 * dx1 * (dx0 + dx1)) * direction[i];
		}
		else if (direction[i] < 0) {
			if (!node->neighbor[i * 2])
				continue;
			double dx;
			if (node->neighbor[i * 2]->isleaf()) {
				x0 = x[node->neighbor[i * 2]->index];
				dx = node->dx + node->neighbor[i * 2]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x0 += x[node->neighbor[i * 2]->child[graph[i * 2][j]]->index];
				}
				x0 /= 4;
				dx = 1.5 * node->dx;
			}
			y[i] = (x1 - x0) / dx * direction[i];
		}
		else if (direction[i] > 0) {
			if (!node->neighbor[i * 2 + 1])
				continue;
			double dx;
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				x2 = x[node->neighbor[i * 2 + 1]->index];
				dx = node->dx + node->neighbor[i * 2 + 1]->dx;
			}
			else {
				for (int j = 0; j < 4; j++) {
					x2 += x[node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index];
				}
				x2 /= 4;
				dx = 1.5 * node->dx;
			}
			y[i] = (x2 - x1) / dx * direction[i];
		}
	}
	return;
}
void PDEOperator::ForwardPNPNS(double* x, double* y) {
	long long int n = mesh->index.size();
	double* v = x;
	double* c0 = v + n;
	double* c1 = c0 + n;
	double* u0 = c1 + n;
	double* u1 = u0 + n;
	double* u2 = u1 + n;
	double* p = u2 + n;
	std::vector<TreeNode*>& nodes = mesh->index;
	for (long long int i = 0; i < n; i++) {
		TreeNode* node = nodes[i];
		NodeType type = node->type;
		double res[3] = { 0, 0, 0 };
		double coeff[6] = { 1, 1, 1, 1, 1, 1 };
		if (type == NodeType::Water) {
			double tmp[3] = { 0 ,0, 0 };
			double tmp1[3] = { 0 ,0, 0 };
			double direction[3] = { u0[i], u1[i], u2[i] };
			this->lapl(node, v, tmp, coeff);
			res[0] = tmp[0] + tmp[1] + tmp[2] + c0[i] * charge_coeff0 + c1[i] * charge_coeff1;
			res[1] = (tmp[0] + tmp[1] + tmp[2]) * c0[i] * z0;
			res[2] = (tmp[0] + tmp[1] + tmp[2]) * c1[i] * z1;

			this->grad(node, v, tmp);
			this->lapl(node, c0, tmp1, coeff);
			res[1] += tmp1[0] + tmp1[1] + tmp1[2];
			this->grad(node, c0, tmp1);
			res[1] += (tmp[0] * tmp1[0] + tmp[1] * tmp1[1] + tmp[2] * tmp1[2]) * z0;
			this->grad_direction(node, c0, tmp1, direction);
			res[1] -= np_ns_coeff0 * (tmp1[0] + tmp1[1] + tmp1[2]);

			this->lapl(node, c1, tmp1, coeff);
			res[2] += tmp1[0] + tmp1[1] + tmp1[2];
			this->grad(node, c1, tmp1);
			res[2] += (tmp[0] * tmp1[0] + tmp[1] * tmp1[1] + tmp[2] * tmp1[2]) * z1;
			this->grad_direction(node, c1, tmp1, direction);
			res[2] -= np_ns_coeff1 * (tmp1[0] + tmp1[1] + tmp1[2]);
			y[i] = res[0];
			y[i + n] = res[1];
			y[i + n * 2] = res[2];
		}
		else if (type == NodeType::ProteinBoundary) {
			double norm[3] = { node->norm[0],  node->norm[1], node->norm[2] };
			double tmp[3] = { 0 ,0, 0 };
			this->grad(node, v, tmp);
			res[1] = (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c0[i] * z0;
			res[2] = (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c1[i] * z1;
			this->grad_direction(node, c0, tmp, norm);
			res[1] += tmp[0] + tmp[1] + tmp[2];
			res[1] -= np_ns_coeff0 * (c0[i] * u0[i] * norm[0] + c0[i] * u1[i] * norm[1] + c0[i] * u2[i] * norm[2]);
			this->grad_direction(node, c1, tmp, norm);
			res[2] += tmp[0] + tmp[1] + tmp[2];
			res[2] -= np_ns_coeff1 * (c1[i] * u0[i] * norm[0] + c1[i] * u1[i] * norm[1] + c1[i] * u2[i] * norm[2]);
			for (int j = 0; j < 6; j++) {
				if (node->neighbor[j]->type == node->type && node->type == NodeType::Protein) {
					coeff[j] = this->diep;
				}
				else if (node->neighbor[j]->type != node->type) {
					coeff[j] = 2 * diep / (1 + diep);
				}
			}
			this->lapl(node, v, tmp, coeff);
			y[i] = tmp[0] + tmp[1] + tmp[2];
			y[i + n] = res[1];
			y[i + n * 2] = res[2];
		}
		else if (type == NodeType::PoreBoundary || type == NodeType::BoundaryN) {
			double norm[3] = { node->norm[0],  node->norm[1], node->norm[2] };
			double tmp[3] = { 0 ,0, 0 };
			this->grad(node, v, tmp);
			res[0] = tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2];
			res[1] = (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c0[i] * z0;
			res[2] = (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c1[i] * z1;
			this->grad(node, c0, tmp);
			res[1] += tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2];
			res[1] -= np_ns_coeff0 * (c0[i] * u0[i] * norm[0] + c0[i] * u1[i] * norm[1] + c0[i] * u2[i] * norm[2]);
			this->grad(node, c1, tmp);
			res[2] += tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2];
			res[2] -= np_ns_coeff1 * (c1[i] * u0[i] * norm[0] + c1[i] * u1[i] * norm[1] + c1[i] * u2[i] * norm[2]);

			y[i] = res[0];
			y[i + n] = res[1];
			y[i + n * 2] = res[2];
		}
		else if (type == NodeType::BoundaryD) {
			y[i] = x[i];
			y[i + n] = c0[i];
			y[i + n * 2] = c1[i];
		}
		else if (type == NodeType::Protein || type == NodeType::WaterVoids) {
			//check the protein or protein boundary or water voids neignbour are same level;
			y[i + n] = c0[i];
			y[i + n * 2] = c1[i];
			for (int j = 0; j < 6; j++) {
				if (node->neighbor[j]->type == node->type && node->type == NodeType::Protein) {
					coeff[j] = this->diep;
				}
				else if (node->neighbor[j]->type != node->type) {
					coeff[j] = 2 * diep / (1 + diep);
				}
			}
			double tmp[3] = { 0 ,0, 0 };
			this->lapl(node, v, tmp, coeff);
			y[i] = tmp[0] + tmp[1] + tmp[2];
		}
	}
}

void PDEOperator::JacobianPNPNS(PNPNS::CSRMatrix& J, double* x) {
	long long int n = mesh->index.size();
	double* v = x;
	double* c0 = v + n;
	double* c1 = c0 + n;
	double* u0 = c1 + n;
	double* u1 = u0 + n;
	double* u2 = u1 + n;
	double* p = u2 + n;
	std::vector<TreeNode*>& nodes = mesh->index;
	std::vector<PNPNS::ETrip> Jlist;
	for (long long int i = 0; i < n; i++) {
		TreeNode* node = nodes[i];
		NodeType type = node->type;
		double res[3] = { 0, 0, 0 };
		double coeff[6] = { 1, 1, 1, 1, 1, 1 };
		if (type == NodeType::Water) {
			double tmp[3] = { 0 ,0, 0 };
			PartialList y;
			this->laplJ(node, y, coeff);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i, it->index, it->value);
				Jlist.emplace_back(i + n, it->index, it->value * c0[i] * z0);
				Jlist.emplace_back(i + 2 * n, it->index, it->value * c1[i] * z1);
				Jlist.emplace_back(i + n, it->index + n, it->value * 1 * c0[it->index]);
				Jlist.emplace_back(i + 2 * n, it->index + 2 * n, it->value * 1 * c1[it->index]);
			}
			Jlist.emplace_back(i, i + n, charge_coeff0 * 1 * c0[i]);
			Jlist.emplace_back(i, i + 2 * n, charge_coeff1 * 1 * c1[i]);
			this->lapl(node, v, tmp, coeff);
			double value = tmp[0] + tmp[1] + tmp[2];
			Jlist.emplace_back(i + n, i + n, value * z0 * 1 * c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, value * z1 * 1 * c1[i]);
			this->grad(node, v, tmp);
			this->gradJ(node, y, tmp);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i + n, it->index + n, it->value * z0 * 1 * c0[it->index]);
				Jlist.emplace_back(i + 2 * n, it->index + 2 * n, it->value * z1 * 1 * c1[it->index]);
			}
			this->grad(node, c0, tmp);
			this->gradJ(node, y, tmp);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i + n, it->index, it->value * z0);
			}
			this->grad(node, c1, tmp);
			this->gradJ(node, y, tmp);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i + 2 * n, it->index, it->value * z1);
			}
			double uvector[3] = {u0[i], u1[i], u2[i]};
			this->gradJ_direction(node, y, uvector);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i + n, it->index + n, -it->value * np_ns_coeff0 * 1 * c0[it->index]);
				Jlist.emplace_back(i + 2 * n, it->index + 2 * n, -it->value * np_ns_coeff1 * 1 * c1[it->index]);
			}
		}
		else if (type == NodeType::ProteinBoundary) {
			double norm[3] = {node->norm[0], node->norm[1], node->norm[2]};
			double tmp[3] = { 0 ,0, 0 };
			for (int j = 0; j < 6; j++) {
				if (node->neighbor[j]->type == node->type && node->type == NodeType::Protein) {
					coeff[j] = this->diep;
				}
				else if (node->neighbor[j]->type != node->type) {
					coeff[j] = 2 * diep / (1 + diep);
				}
			}
			PartialList y;
			this->laplJ(node, y, coeff);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i, it->index, it->value);
			}
			Jlist.emplace_back(i, i + n, charge_coeff0 * 1 * c0[i]);
			Jlist.emplace_back(i, i + 2 * n, charge_coeff1 * 1 * c1[i]);
			this->gradJ2(node, y, norm);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i + n, it->index + n, it->value * 1 * c0[it->index]);
				Jlist.emplace_back(i + 2 * n, it->index + 2 * n, it->value * 1 * c1[it->index]);
			}
			this->gradJ2(node, y, norm);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i + n, it->index, it->value * z0 * c0[i]);
				Jlist.emplace_back(i + 2 * n, it->index, it->value * z1 * c1[i]);
			}
			this->grad(node, v, tmp);
			Jlist.emplace_back(i + n, i + n, z0 * (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, z1 * (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c1[i]);
			Jlist.emplace_back(i + n, i + n, -np_ns_coeff0 * (u0[i] * norm[0] + u1[i] * norm[1] + u2[i] * norm[2]) * c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, -np_ns_coeff1 * (u0[i] * norm[0] + u1[i] * norm[1] + u2[i] * norm[2]) * c1[i]);
		}
		else if (type == NodeType::PoreBoundary || type == NodeType::BoundaryN) {
			double norm[3] = { node->norm[0], node->norm[1], node->norm[2] };
			double tmp[3] = { 1 , 1, 1 };
			PartialList y;
			this->gradJ(node, y, norm);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i, it->index, it->value);
				Jlist.emplace_back(i + n, it->index, it->value * z0 * c0[i]);
				Jlist.emplace_back(i + 2 * n, it->index, it->value * z1 * c1[i]);
				Jlist.emplace_back(i + n, it->index + n, it->value * c0[it->index]);
				Jlist.emplace_back(i + 2 * n, it->index + 2 * n, it->value * c1[it->index]);
			}
			this->grad(node, v, tmp);
			Jlist.emplace_back(i + n, i + n, z0 * (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, z1 * (tmp[0] * norm[0] + tmp[1] * norm[1] + tmp[2] * norm[2]) * c1[i]);
			Jlist.emplace_back(i + n, i + n, -np_ns_coeff0 * (u0[i] * norm[0] + u1[i] * norm[1] + u2[i] * norm[2]) * c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, -np_ns_coeff1 * (u0[i] * norm[0] + u1[i] * norm[1] + u2[i] * norm[2]) * c1[i]);
		}
		else if (type == NodeType::BoundaryD) {
			Jlist.emplace_back(i + n, i + n, c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, c1[i]);
			Jlist.emplace_back(i, i, 1);
		}
		else if (type == NodeType::Protein || type == NodeType::WaterVoids) {
			Jlist.emplace_back(i + n, i + n, c0[i]);
			Jlist.emplace_back(i + 2 * n, i + 2 * n, c1[i]);
			for (int j = 0; j < 6; j++) {
				if (node->neighbor[j]->type == node->type && node->type == NodeType::Protein) {
					coeff[j] = diep;
				}
				else if (node->neighbor[j]->type != node->type) {
					coeff[j] = 2 * diep / (1 + diep);
				}
			}
			PartialList y;
			this->laplJ(node, y, coeff);
			for (auto it = y.begin(); it != y.end(); ++it) {
				Jlist.emplace_back(i, it->index, it->value);
			}
		}
	}
	for (long long int i = 0; i < Jlist.size(); i++) {
		if (Jlist[i].col() < 0 || Jlist[i].row() < 0) {
			long long int n = Jlist[i].row();
			if (n < mesh->index.size()) {
				long long int ntype = mesh->index[n]->type;
				std::cout << "row: " << Jlist[i].row() << " col: " << Jlist[i].col() << " type: "<< ntype << " index: "<< mesh->index[n]->index << std::endl;
			}
		}
	}
	J.setFromTriplets(Jlist.begin(), Jlist.end());
}
void PDEOperator::gradJ(TreeNode* node, PartialList& y, double* coeffs) {
	y.clear();
	long long int id0 = node->index;
	double value0 = 0;
	for (int i = 0; i < 3; i++) {
		if (node->neighbor[i * 2] && node->neighbor[i * 2 + 1]) {
			if (node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx0 = 2 * node->dx;
				//double dx1 = 2 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
			}
			else if(!node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = 1.5 * node->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx1 = 2 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else if (node->neighbor[i * 2]->isleaf() && !node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				//double dx0 = 2 * node->dx;
				double dx1 = 1.5 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
			}
			else {
				double coeff = 1.0 / 3 / node->dx * coeffs[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, -coeff / 4);
				}
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (node->neighbor[i * 2]) {
			if (node->neighbor[i * 2]->isleaf()) {
				double coeff = 1.0 / (node->dx + node->neighbor[i * 2]->dx) * coeffs[i];
				//double coeff = 1.0 / (2 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
			}
			else {
				double coeff = 1.0 / (1.5 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (node->neighbor[i * 2 + 1]) {
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				double coeff = -1.0 / (node->dx + node->neighbor[i * 2 + 1]->dx) * coeffs[i];
				//double coeff = -1.0 / (2 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else {
				double coeff = -1.0 / (1.5 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
	}
	if (value0 != 0)
		y.emplace_front(id0, value0);
	return;
};
void PDEOperator::gradJ2(TreeNode* node, PartialList& y, double* coeffs) {
	y.clear();
	long long int id0 = node->index;
	double value0 = 0;
	for (int i = 0; i < 3; i++) {
		if (node->neighbor[i * 2]->type == NodeType::Water && node->neighbor[i * 2 + 1]->type == NodeType::Water) {
			if (node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx0 = 2 * node->dx;
				//double dx1 = 2 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
			}
			else if (!node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = 1.5 * node->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx1 = 2 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else if (node->neighbor[i * 2]->isleaf() && !node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				//double dx0 = 2 * node->dx;
				double dx1 = 1.5 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * coeffs[i];
			}
			else {
				double coeff = 1.0 / 3 / node->dx * coeffs[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, -coeff / 4);
				}
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (node->neighbor[i * 2]->type == NodeType::Water) {
			if (node->neighbor[i * 2]->isleaf()) {
				double coeff = 1.0 / (node->dx + node->neighbor[i * 2]->dx) * coeffs[i];
				//double coeff = 1.0 / (2 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
			}
			else {
				double coeff = 1.0 / (1.5 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (node->neighbor[i * 2 + 1]->type == NodeType::Water) {
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				double coeff = -1.0 / (node->dx + node->neighbor[i * 2 + 1]->dx) * coeffs[i];
				//double coeff = -1.0 / (2 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else {
				double coeff = -1.0 / (1.5 * node->dx) * coeffs[i];
				value0 += coeff;
				coeff = -coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
	}
	if (value0 != 0)
		y.emplace_front(id0, value0);
	return;
};
void PDEOperator::laplJ(TreeNode* node, PartialList& y, double* coeffs) {
	y.clear();
	long long int id0 = node->index;
	double value0 = 0;
	for (int i = 0; i < 3; i++) {
		if (node->neighbor[i * 2] && node->neighbor[i * 2 + 1]) {
			if (node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx0 = 2 * node->dx;
				//double dx1 = 2 * node->dx;
				double coeff = 2.0 / (dx0 * (dx0 + dx1)) * coeffs[i * 2];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = 2.0 / (dx1 * (dx0 + dx1)) * coeffs[i * 2 + 1];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
				value0 -= (coeffs[i * 2] + coeffs[i * 2 + 1]) / dx0 / dx1;
			}
			else if (!node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = 1.5 * node->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx1 = 2 * node->dx;
				double coeff = 2.0 / (dx0 * (dx0 + dx1)) * coeffs[2 * i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				coeff = 2.0 / (dx1 * (dx0 + dx1)) * coeffs[2 * i + 1];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
				value0 -= (coeffs[i * 2] + coeffs[i * 2 + 1]) / (dx1 * dx0);
			}
			else if (node->neighbor[i * 2]->isleaf() && !node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				//double dx0 = 2 * node->dx;
				double dx1 = 1.5 * node->dx;
				double coeff = 2.0 / (dx0 * (dx0 + dx1)) * coeffs[2 * i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = 2.0 / (dx1 * (dx0 + dx1)) * coeffs[2 * i + 1];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 -= (coeffs[i * 2] + coeffs[i * 2 + 1]) / (dx1 * dx0);
			}
			else {
				double dx = node->dx * 1.5;
				double coeff = 1.0 / dx / dx;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4 * coeffs[2 * i]);
				}
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4 * coeffs[2 * i + 1]);
				}
				value0 -= (coeffs[i * 2] + coeffs[i * 2 + 1]) * coeff;
			}
		}
		else if (node->neighbor[i * 2]) {
			if (node->neighbor[i * 2]->isleaf()) {
				double dx = (node->dx + node->neighbor[i * 2]->dx);
				//double dx = 2 * node->dx;
				double coeff = 1.0 / dx / dx * coeffs[2 * i];
				value0 -= coeff;
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
			}
			else {
				double dx = 1.5 * node->dx;
				double coeff = 1.0 / dx / dx * coeffs[2 * i];
				value0 -= coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (node->neighbor[i * 2 + 1]) {
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				double dx = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx = 2 * node->dx;
				double coeff = 1.0 / dx / dx * coeffs[2 * i + 1];
				value0 -= coeff;
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else {
				double dx = 1.5 * node->dx;
				double coeff = 1.0 / dx / dx * coeffs[2 * i + 1];
				value0 -= coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
	}
	if (value0 != 0)
		y.emplace_front(id0, value0);
	
	return;
};
void PDEOperator::gradJ_direction(TreeNode* node, PartialList& y, double* direction) {
	y.clear();
	long long int id0 = node->index;
	double value0 = 0;
	for (int i = 0; i < 3; i++) {
		if (direction[i] == 0) {
			if (node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				//double dx0 = 2 * node->dx;
				//double dx1 = 2 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
			}
			else if (!node->neighbor[i * 2]->isleaf() && node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = 1.5 * node->dx;
				//double dx1 = 2 * node->dx;
				double dx1 = node->dx + node->neighbor[i * 2 + 1]->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else if (node->neighbor[i * 2]->isleaf() && !node->neighbor[i * 2 + 1]->isleaf()) {
				double dx0 = node->dx + node->neighbor[i * 2]->dx;
				//double dx0 = 2 * node->dx;
				double dx1 = 1.5 * node->dx;
				double coeff = -dx1 * dx1 / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
				coeff = dx0 * dx0 / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
				value0 += (dx1 * dx1 - dx0 * dx0) / (dx1 * dx0 * (dx0 + dx1)) * direction[i];
			}
			else {
				double coeff = 1.0 / (3 * node->dx) * direction[i];
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, -coeff / 4);
				}
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (direction[i] < 0) {
			if (node->neighbor[i * 2]->isleaf()) {
				double coeff = 1.0 / (node->dx + node->neighbor[i * 2]->dx) * direction[i];
				//double coeff = 1.0 / (node->dx * 2) * direction[i];
				value0 += coeff;
				coeff = -coeff;
				y.emplace_back(node->neighbor[i * 2]->index, coeff);
			}
			else {
				double coeff = 1.0 / (1.5 * node->dx) * direction[i];
				value0 += coeff;
				coeff = -coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2]->child[graph[i * 2][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
		else if (direction[i] > 0) {
			if (node->neighbor[i * 2 + 1]->isleaf()) {
				double coeff = -1.0 / (node->dx + node->neighbor[i * 2 + 1]->dx) * direction[i];
				//double coeff = -1.0 / (2 * node->dx) * direction[i];
				value0 += coeff;
				coeff = -coeff;
				y.emplace_back(node->neighbor[i * 2 + 1]->index, coeff);
			}
			else {
				double coeff = -1.0 / (1.5 * node->dx) * direction[i];
				value0 += coeff;
				coeff = -coeff;
				for (int j = 0; j < 4; j++) {
					long long int id = node->neighbor[i * 2 + 1]->child[graph[i * 2 + 1][j]]->index;
					y.emplace_back(id, coeff / 4);
				}
			}
		}
	}
	if (value0 != 0)
		y.emplace_front(id0, value0);
	
	return;
};

void PDEOperator::righthand(double* x, double v0, double v1, double surfacedensity) {
	double charge_coeff_volume = -ekbt * ee * 1e10 / eps0;
	double charge_coeff_surface = -ekbt * 1e-10 / eps0;
	long long int n = mesh->index.size();
	std::vector<TreeNode*>& nodes = mesh->index;
	for (long long int i = 0; i < n; i++) {
		TreeNode* node = nodes[i];
		if (node->isleaf() && node->type == NodeType::PoreBoundary) {
			x[i] -= surfacedensity * charge_coeff_surface;
		}
		else if (node->isleaf() && node->type == NodeType::BoundaryD) {
			if (node->z > 0) {
				x[i] -= v1 * ekbt;
			}
			else {
				x[i] -= v0 * ekbt;
			}
			x[i + n] -= 1;
			x[i + 2 * n] -= 1;
		}
		else if (node->isleaf()) {
			x[i] -= node->chargeDensity * charge_coeff_volume;
		}
	}

};