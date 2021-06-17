#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <pwd.h>
#include <cmath>
#include <boost/property_tree/ptree.hpp>

using namespace std;

vector<vector<vector<double>>> cut_isochrone_into_branches(vector<vector<double>> &isochrone){
    vector<vector<vector<double>>> branches;
    vector<vector<double>> branch;
    for(auto const &elem: isochrone){
        if(elem[0] < 1){
            branch.push_back(elem);
        }
        else{
            branch.push_back(elem);
            branches.push_back(branch);
            branch.clear();
        }
    }
    return branches;
}

double max(double x, double y){
    if(x > y){
        return x;
    }
    else{
        return y;
    }
}

double min(double x, double y){
    if(x < y){
        return x;
    }
    else{
        return y;
    }
}

bool intersect(double ax1, double ay1, double ax2, double ay2, double bx1, double by1, double bx2, double by2) {
    // segment_a = [[ax1, ay1], [ax2, ay2]], segment_b = [[bx1, by1], [bx2, by2]]
    double Aa, Ab, Ba, Bb, x_intersec;
    // Check for overlap: If the largest x coordinate of segment a is smaller than the smallest x coordinate
    // of segment b then there can be no intersection. (same for y)
    if (max(ax1, ax2) < min(bx1, bx2)){
        return false;
    }
    if (max(ay1, ay2) < min(by1, by2)){
        return false;
    }

    // If the is a mutual interval calculate the x coordinate of that intersection point and check if it is in the interval.
    // Calculate fa(a) = Aa*x + Ba = y and fb(x) = Ab*x + Bb = y
    Aa = (ay1 - ay2) / (ax1 - ax2);  // slope of segment a
    Ab = (by1 - by2) / (bx1 - bx2);  // slope of segment b
    Ba = ay1 - Aa * ax1;  // y intercept of segment a
    Bb = by1 - Ab * bx1;  // y intercep of segment b
    x_intersec = (Bb - Ba) / (Aa - Ab);  // x coordinate of intersection point
    if (x_intersec < max(min(ax1, ax2), min(bx1, bx2)) or x_intersec > min(max(ax1, ax2), max(bx1, bx2))){
        return false;
    }
    else{
        return true;
    }
}

int main(int argc, char *argv[]) {

    const int run = atoi(argv[1]);
    float mu = 5.0;
    float D = 1.0;
    float tau_a = 2.0;
    float delta_a = 1.0;
    double pi = 3.14159265;
    double phase = 2.*pi/2.;
    double dt = 0.0001;

    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    std::string path;
    path = "../out/";
    char parameters[200];
    std::sprintf(parameters, "mu%.1f_taua%.1f_delta%.1f_D%.2f_phase%.2f_run%d", mu, tau_a, delta_a, D, phase, run);
    std::string out_file_isi;
    std::string out_file_ipi;
    out_file_isi = path + "ISIs_alif_" + parameters + ".dat";
    out_file_ipi = path + "IPIs_alif_" + parameters + ".dat";
    std::ofstream file_isi;
    std::ofstream file_ipi;
    file_isi.open(out_file_isi);
    file_ipi.open(out_file_ipi);
    if (!file_isi.is_open()) {
        std::cout << "Could not open file at: " << out_file_isi << std::endl;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }
    if (!file_ipi.is_open()) {
        std::cout << "Could not open file at: " << out_file_ipi << std::endl;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }


    double v0 = 0;
    double a0 = 2;


    // The following two lines and the latter two are actually redundant but illustrate what is happening
    int max_spikes = 2000;
    std::vector<double> ISIs(max_spikes);
    std::vector<double> IPIs(max_spikes);
    std::vector<double> v_IPI_pass;
    std::vector<double> a_IPI_pass;

    // Calculate interspike intervals
    bool fired = false;
    int spike_count = 0;
    double t_since_spike = 0;

    double v = 0;
    double a = 0;
    double xi = 0;
    const double mean = 0.0;
    const double stddev = 1.0;
    std::random_device rd;
    std::mt19937 generator(rd());
    //Better seed from random_device instead of clock in case one runs many simulations in a short periode of time
    std::normal_distribution<double> dist(mean, stddev);

    while(spike_count < 100){
        xi = dist(generator);
        v += (mu  - v - a) * dt + sqrt(2 * D * dt)*xi;
        a += (-a/tau_a) * dt;
        if (v > 1.0){
            v = 0;
            a += delta_a;
            spike_count += 1;
        }
    }

    while(spike_count < max_spikes){
        xi = dist(generator);
        v += (mu  - v - a) * dt + sqrt(2 * D * dt)*xi;
        a += (-a/tau_a) * dt;
        t_since_spike += dt;
        if (v > 1.0){
            v = 0;
            a += delta_a;
            ISIs[spike_count] = t_since_spike;
            spike_count += 1;
            t_since_spike = 0;
        }
    }

    // Print interspike intervals to file
    for (auto isi: ISIs) {
        file_isi << isi << "\n";
    }

    // Load isochrone
    std::vector<std::vector<double>> isochrone;
    char phi[200];
    std::sprintf(phi, "%.2f", phase);
    std::string in_file_isochrone = std::string(homedir) + "/Data/isochrones/isochrone_mu5.0_taua2.0_delta1.0_phase" + phi + ".dat";
    std::ifstream file_isochrone(in_file_isochrone, std::ios_base::in);

    std::string line;
    while (std::getline(file_isochrone, line)){
        std::istringstream iss(line);
        double viso, aiso;
        if(!(iss >> viso >> aiso)) {
            break;
        }
        else{
            isochrone.push_back({viso, aiso});
        }
    }
    // Cut isochrone into branches
    vector<vector<vector<double>>> branches = cut_isochrone_into_branches(isochrone);

    // Calculat inter phase intervals
    spike_count = 0;
    t_since_spike = 0;
    int number_branch_to_be_passed = 1;

    bool passed_branch = false;
    vector<vector<double>> branch_to_be_passed;
    double v_before, a_before, v_after, a_after;

    vector<double> p_before, p_after;
    double v_iso_before, a_iso_before, v_iso_after, a_iso_after;

    while(spike_count < max_spikes){
        v_before = v;
        a_before = a;

        xi = dist(generator);
        v += (mu  - v - a) * dt + sqrt(2 * D * dt)*xi;
        a += (-a/tau_a) * dt;
        if (v > 1.0){
            v = 0;
            a += delta_a;
            fired = true;
        }
        else{
            fired = false;
        }
        //fired = alif.integrate_stochastic();

        if(not fired){
            v_after = v;
            a_after = a;
        }
        else{
            v_after = 1.;
            a_after = a - delta_a;
            number_branch_to_be_passed += 1;
        }

        t_since_spike += dt;

        if(number_branch_to_be_passed < 0){
            branch_to_be_passed = {{-1, 0}, {1, 0}};
        }
        else{
            branch_to_be_passed = branches[number_branch_to_be_passed];
        }

        for(int i=0; i < branch_to_be_passed.size()-1; i++){
            p_before = branch_to_be_passed[i];
            p_after = branch_to_be_passed[i+1];
            v_iso_before = p_before[0];
            a_iso_before = p_before[1];
            v_iso_after = p_after[0];
            a_iso_after = p_after[1];
            if (std::abs(v_iso_after - v_iso_before) > 0.5){
                continue;
            }
            passed_branch = intersect(v_before, a_before, v_after, a_after, v_iso_before, a_iso_before, v_iso_after, a_iso_after);
            if(passed_branch){
                cout << spike_count << " " << number_branch_to_be_passed << endl;
                IPIs[spike_count] = t_since_spike;
                spike_count += 1;
                number_branch_to_be_passed -= 1;
                v_IPI_pass.push_back(v);
                a_IPI_pass.push_back(a);
                t_since_spike = 0;
            }
        }
    }

    for(int i = 0; i < IPIs.size(); i++){
        double ipi = IPIs[i];
        double v_print = v_IPI_pass[i];
        double a_print = a_IPI_pass[i];
        file_ipi << ipi << " " << v_print << " " << a_print << " " << "\n";
    }
}
