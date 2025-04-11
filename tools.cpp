#include "tools.h"
#include <fstream>

void printA(PNPNS::CSRMatrix& A, int i) {
    for (int col = 0; col < A.cols(); ++col) {
        for (int row = 0; row < A.rows(); row++) {
            std::cout << "[" << col << "] ";
        }
        std::cout << std::endl;
    }
}

bool checkCSRMatrix(PNPNS::CSRMatrix& A) {
    for (int i = 0; i < A.cols(); i++) {
        if (A.coeff(i, i) == 0) {
            std::cerr << "main diagonal at (" << i << ") 0";
            //return false;
        }
        else {
            //  std::cout << "A: " << A.coeff(i, i) << "at: " << i << std::endl;
        }
    }
    return true;
}

double donutdistance(double r, double o, double off, double a, double b, double c) {
    double d = pow((off - sqrt(a * a + b * b)), 2) + (c - o) * (c - o);
    return sqrt(d);
}

bool cubicwithdonut(double r, double o, double off, double x, double y, double z, double dx) {
    double a = x - dx;
    double b = y - dx;
    double c = z - dx;
    double d0 = donutdistance(r, o, off, a, b, c);
    a = x - dx;
    b = y + dx;
    c = z - dx;
    double d1 = donutdistance(r, o, off, a, b, c);
    a = x + dx;
    b = y - dx;
    c = z - dx;
    double d2 = donutdistance(r, o, off, a, b, c);
    a = x + dx;
    b = y + dx;
    c = z - dx;
    double d3 = donutdistance(r, o, off, a, b, c);
    a = x - dx;
    b = y - dx;
    c = z + dx;
    double d4 = donutdistance(r, o, off, a, b, c);
    a = x - dx;
    b = y + dx;
    c = z + dx;
    double d5 = donutdistance(r, o, off, a, b, c);
    a = x + dx;
    b = y - dx;
    c = z + dx;
    double d6 = donutdistance(r, o, off, a, b, c);
    a = x + dx;
    b = y + dx;
    c = z + dx;
    double d7 = donutdistance(r, o, off, a, b, c);
    if (d0 > r && d1 > r && d2 > r && d3 > r && d4 > r && d5 > r && d6 > r && d7 > r) {
        return false;
    }
    else if (d0 < r && d1 < r && d2 < r && d3 < r && d4 < r && d5 < r && d6 < r && d7 < r) {
        return false;
    }
    return true;
}


void donutnormal(double r, double o, double off, double a, double b, double c, double* res) {
    double d = sqrt(a * a + b * b);
    res[0] = a - a * off / d;
    res[1] = b - b * off / d;
    res[2] = c - o;
    d = sqrt(res[0] * res[0] + res[1] * res[1] + res[2] * res[2]);
    res[0] /= d;
    res[1] /= d;
    res[2] /= d;
    return;
}
int pointnanopore(double x, double y, double z, double r, double l) {
    double r0 = r * 0.2;
    double d = sqrt(x * x + y * y);
    if (z > l || d < r)
        return 1;
    if (z == l && d >= (r + r0) || d == r && z <= l - r0)
        return 0;
    if (z < l && d >= (r + r0) || d < (r + r0) && d > r && z < l - r0) {
        return -1;
    }
    d = donutdistance(r0, l - r0, r + r0, x, y, z);
    if (d > r0)
        return 1;
    if (d < r0)
        return -1;
    return 0;
}
int pointnanopore2(double x, double y, double z, double r, double l, double debye) {
    double dr = sqrt(x * x + y * y);
    double d = 0;
    int res = 0;
    if (dr >= 1.2 * r) {
        d = z - l;
    }
    else if (z <= l - 0.2 * r) {
        d = r - dr;
    }
    else {
        d = donutdistance(r * 0.2, l - 0.2 * r, r * 1.2, x, y, z);
    }
    if (d <= debye && d >= -debye) {
        return 0;
    }
    else if (d < 0) {
        return -1;
    }
    else
        return 1;
}
int cubicwithpore(double x, double y, double z, double dx, double r, double l, double debye) {
    double x0; double y0; double z0;
    if (z - dx >= 0 || z + dx <= 0) {
        x0 = x - dx; y0 = y - dx; z0 = abs(z - dx);
        int c0 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c0 == 0)
            return 0;
        x0 = x - dx; y0 = y + dx; z0 = abs(z - dx);
        int c1 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c1 == 0)
            return 0;
        x0 = x + dx; y0 = y - dx; z0 = abs(z - dx);
        int c2 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c2 == 0)
            return 0;
        x0 = x + dx; y0 = y + dx; z0 = abs(z - dx);
        int c3 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c3 == 0)
            return 0;
        x0 = x - dx; y0 = y - dx; z0 = abs(z + dx);
        int c4 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c4 == 0)
            return 0;
        x0 = x - dx; y0 = y + dx; z0 = abs(z + dx);
        int c5 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c5 == 0)
            return 0;
        x0 = x + dx; y0 = y - dx; z0 = abs(z + dx);
        int c6 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c6 == 0)
            return 0;
        x0 = x + dx; y0 = y + dx; z0 = abs(z + dx);
        int c7 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c7 == 0)
            return 0;
        if (c0 == 1 && c1 == 1 && c2 == 1 && c3 == 1 &&
            c4 == 1 && c5 == 1 && c6 == 1 && c7 == 1)
            return 1; 
        else if (c0 == -1 && c1 == -1 && c2 == -1 && c3 == -1 &&
            c4 == -1 && c5 == -1 && c6 == -1 && c7 == -1) {
            return -1;
        }
        else {
            return 0;
        }
    }
    else {
        int res0;
        x0 = x - dx; y0 = y - dx; z0 = abs(z + dx);
        int c0 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c0 == 0)
            return 0;
        x0 = x - dx; y0 = y + dx; z0 = abs(z + dx);
        int c1 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c1 == 0)
            return 0;
        x0 = x + dx; y0 = y - dx; z0 = abs(z + dx);
        int c2 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c2 == 0)
            return 0;
        x0 = x + dx; y0 = y + dx; z0 = abs(z + dx);
        int c3 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c3 == 0)
            return 0;
        x0 = x - dx; y0 = y - dx; z0 = 0;
        int c4 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c4 == 0)
            return 0;
        x0 = x - dx; y0 = y + dx; z0 = 0;
        int c5 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c5 == 0)
            return 0;
        x0 = x + dx; y0 = y - dx; z0 = 0;
        int c6 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c6 == 0)
            return 0;
        x0 = x + dx; y0 = y + dx; z0 = 0;
        int c7 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c7 == 0)
            return 0;
        if (c0 == 1 && c1 == 1 && c2 == 1 && c3 == 1 &&
            c4 == 1 && c5 == 1 && c6 == 1 && c7 == 1)
            res0 = 1;
        else if (c0 == -1 && c1 == -1 && c2 == -1 && c3 == -1 &&
            c4 == -1 && c5 == -1 && c6 == -1 && c7 == -1) {
            res0 = -1;
        }
        else {
            res0 = 0;
        }
        int res1;
        x0 = x - dx; y0 = y - dx; z0 = 0;
        c0 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c0 == 0)
            return 0;
        x0 = x - dx; y0 = y + dx; z0 = 0;
        c1 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c1 == 0)
            return 0;
        x0 = x + dx; y0 = y - dx; z0 = 0;
        c2 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c2 == 0)
            return 0;
        x0 = x + dx; y0 = y + dx; z0 = 0;
        c3 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c3 == 0)
            return 0;
        x0 = x - dx; y0 = y - dx; z0 = abs(z + dx);
        c4 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c4 == 0)
            return 0;
        x0 = x - dx; y0 = y + dx; z0 = abs(z + dx);
        c5 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c5 == 0)
            return 0;
        x0 = x + dx; y0 = y - dx; z0 = abs(z + dx);
        c6 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c6 == 0)
            return 0;
        x0 = x + dx; y0 = y + dx; z0 = abs(z + dx);
        c7 = pointnanopore2(x0, y0, z0, r, l / 2, debye);
        if (c7 == 0)
            return 0;
        if (c0 == 1 && c1 == 1 && c2 == 1 && c3 == 1 &&
            c4 == 1 && c5 == 1 && c6 == 1 && c7 == 1)
            res1 = 1;
        else if (c0 == -1 && c1 == -1 && c2 == -1 && c3 == -1 &&
            c4 == -1 && c5 == -1 && c6 == -1 && c7 == -1) {
            res1 = -1;
        }
        else {
            res1 = 0;
        }
        if (res0 == 1 && res1 == 1)
            return 1;
        else if (res0 == -1 && res1 == -1)
            return -1;
        else
            return 0;
    }
}

void readPQR(PNPNS::Atoms& atoms, std::string& fn, double probe, double mx, double my, double mz) {
    //read pqrt
    std::ifstream file(fn);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << fn << std::endl;
        return;
    }
    std::string line;
    PNPNS::double3 pos_min;
    PNPNS::double3 pos_max;
    PNPNS::double3 pos_mean = { 0, 0, 0 };
    int n_atom = 0;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM") {
            double x = std::stof(line.substr(30, 8));
            atoms.x.push_back(x);
            x = std::stof(line.substr(38, 8));
            atoms.y.push_back(x);
            x = std::stof(line.substr(46, 8));
            atoms.z.push_back(x);
            atoms.charge.push_back(std::stof(line.substr(54, 7)));
            x = std::stof(line.substr(61, 7));
            x = x + probe;
            atoms.radius.push_back(x);
            pos_mean.x += atoms.x.back();
            pos_mean.y += atoms.y.back();
            pos_mean.z += atoms.z.back();
            double max_radius = atoms.radius.back();
            if (n_atom == 0) {
                pos_min = { atoms.x.back() - max_radius, atoms.y.back() - max_radius, atoms.z.back() - max_radius };
                pos_max = { atoms.x.back() + max_radius, atoms.y.back() + max_radius, atoms.z.back() + max_radius };
            }
            else {
                pos_min.x = (atoms.x.back() - max_radius < pos_min.x) ? atoms.x.back() - max_radius : pos_min.x;
                pos_min.y = (atoms.y.back() - max_radius < pos_min.y) ? atoms.y.back() - max_radius : pos_min.y;
                pos_min.z = (atoms.z.back() - max_radius < pos_min.z) ? atoms.z.back() - max_radius : pos_min.z;
                pos_max.x = (atoms.x.back() + max_radius > pos_max.x) ? atoms.x.back() + max_radius : pos_max.x;
                pos_max.y = (atoms.y.back() + max_radius > pos_max.y) ? atoms.y.back() + max_radius : pos_max.y;
                pos_max.z = (atoms.z.back() + max_radius > pos_max.z) ? atoms.z.back() + max_radius : pos_max.z;
            }
            n_atom++;
        }
    }
    file.close();
    pos_mean.x /= n_atom;
    pos_mean.y /= n_atom;
    pos_mean.z /= n_atom;
    pos_min.x = pos_min.x - pos_mean.x + mx;
    pos_min.y = pos_min.y - pos_mean.y + my;
    pos_min.z = pos_min.z - pos_mean.z + mz;
    pos_max.x = pos_max.x - pos_mean.x + mx;
    pos_max.z = pos_max.y - pos_mean.y + my;
    pos_max.z = pos_max.z - pos_mean.z + mz;
    for (int i = 0; i < n_atom; i++) {
        atoms.x[i] = atoms.x[i] - pos_mean.x + mx;
        atoms.y[i] = atoms.y[i] - pos_mean.y + my;
        atoms.z[i] = atoms.z[i] - pos_mean.z + mz;
    }
    std::cout << "protein coordinate range: [" << pos_min.x << ", " << pos_min.y << ", " << pos_min.z << "] ["
        << pos_max.x << ", " << pos_max.y << ", " << pos_max.z << "]" << std::endl;
}