#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>

// Custom data structure to make storing double precision 2D vector data easier
struct double2{
    double x, y;
};

// Custom data structure to make storing integer 2D vector data easier
struct int2{
    int x, y;
};

// Custom data structure to make storing single precision 2D vector data easier
struct float2{
    float x, y;
};

// Helper function to output the "galaxies" and randoms to files
void writeFile(std::string fileName, std::vector<double2> points) {
    std::ofstream fout(fileName);
    fout.precision(15);
    for (int i = 0; i < points.size(); ++i) {
        fout << points[i].x << " " << points[i].y << "\n";
    }
    fout.close();
}

// Helper function to output the grid values to files
void writeFile(std::string fileName, std::vector<int> grid, int2 N) {
    std::ofstream fout(fileName);
    for (int i = 0; i < N.x; ++i) {
        for (int j = 0; j < N.y; ++j) {
            int index = j + N.y*i;
            fout << i << " " << j << " " << grid[index] << "\n";
        }
    }
    fout.close();
}

// Main function, all C++ programs start here
int main() {
    // Set up the grid dimensions and sizes for this simple example
    int2 N = {8,8};
    double2 L = {8,8};
    double2 dr = {L.x/N.x, L.y/N.y};
    
    // This is just setting up a random number generator in C++ with two different distributions that can 
    // be drawn from; a normal distribution to mimic galaxy clusters, and a uniform distribution for the randoms
    std::random_device seeder; // System entropy random number generator
    std::mt19937_64 gen(seeder()); // Mersenne twister random number generator seeded with system entropy
    std::normal_distribution<double> nDist(0.0,0.5);
    std::uniform_real_distribution<double> uDist(0.0,1.0);
    
    // Numbers of "galaxies" and randoms to generate
    int N_gals = 100;
    int N_rans = 1000;
    
    // Arrays to store the "galaxies" and randoms (i.e. their positions)
    std::vector<double2> gals;
    std::vector<double2> rans;
    
    // This just generates the "galaxies" by first setting the center of the cluster, and then drawing shifts
    // away from that center from a Gaussian distribution. The if statements implement periodic boundary 
    // conditions so that all galaxies have positions between 0 and 8 in x and y.
    for (int clump = 0; clump < 4; ++clump) {
        double2 center = {uDist(gen)*L.x, uDist(gen)*L.y};
        for (int i = 0; i < N_gals/4; ++i) {
            double2 gal = {center.x + nDist(gen), center.y + nDist(gen)};
            if (gal.x >= L.x) gal.x -= L.x;
            if (gal.x < 0) gal.x += L.x;
            if (gal.y >= L.y) gal.y -= L.y;
            if (gal.y < 0) gal.y += L.y;
            gals.push_back(gal);
        }
    }
    
    // Generates the uniform random points
    for (int i = 0; i < N_rans; ++i) {
        double2 ran = {uDist(gen)*L.x, uDist(gen)*L.y};
        rans.push_back(ran);
    }
    
    // Output the lists to files.
    writeFile("gals.dat", gals);
    writeFile("rans.dat", rans);
    
    // Create the grid arrays to bin the "galaxies and randoms to
    std::vector<int> d_gals(N.x*N.y);
    std::vector<int> d_rans(N.x*N.y);
    
    // Bin the galaxies
    for (int i = 0; i < gals.size(); ++i) {
        int x = int(gals[i].x/dr.x);
        int y = int(gals[i].y/dr.y);
        // This index variable is because in C/C++ we use flattened arrays when doing FFTs, in Python, you
        // should be able to simply have d_gal[x][y] += 1
        int index = y + N.y*x;
        
        d_gals[index] += 1;
    }
    
    // Bin the randoms
    for (int i = 0; i < rans.size(); ++i) {
        int x = int(rans[i].x/dr.x);
        int y = int(rans[i].y/dr.y);
        int index = y + N.y*x;
        
        d_rans[index] += 1;
    }
    
    // Write the grid arrays to files
    writeFile("galsGrid.dat", d_gals, N);
    writeFile("ransGrid.dat", d_rans, N);
    
    // End program
    return 0;
}
        
