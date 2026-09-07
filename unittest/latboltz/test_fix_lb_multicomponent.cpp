/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 Contributing authors: Abhin Suresh (University of Delaware)
------------------------------------------------------------------------- */

// unit tests for fix style multicomponent fluid dynamics

#include "test_main.h"
#include "yaml_writer.h"

#include "gtest/gtest.h"

#include "input.h"
#include "fix.h"
#include "fix_lb_multicomponent.h"
#include "lammps.h"
#include "modify.h"

#include <iostream>
#include <mpi.h>
#include <string>

using namespace LAMMPS_NS;

LAMMPS *init_lammps()
{
    const char *args[] = {"FixLBMulticomponent", "-log", "none", 
                          "-echo", "screen", "-nocite"};
    
    char **argv = (char **)args;
    int argc = sizeof(args) / sizeof(char *);

    LAMMPS *lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);

    // utility lambda to improve readability
    auto command = [&](const std::string &line){
        lmp->input->one(line); 
    };
    
    command("region fluid block -16 16 -16 16 -3 3");
    command("create_box 0 fluid");
    command("timestep 1.0");
    command("fix mcmp all lb/multicomponent "
            "1 0.166667 1.0 D3Q19 dx 1 "
            "C1 0.333333 C2 0.333333 C3 0.333334 "
            "kappa1 0.01 kappa2 0.02 kappa3 0.05 "
            "init mixture");

    return lmp;
}

void run_lammps(LAMMPS *lmp)
{
    lmp->input->one("run 100");
}

void generate_yaml_file(const char *outfile, const TestConfig &config)
{
    LAMMPS *lmp = init_lammps();

    YamlWriter writer(outfile);

    write_yaml_header(&writer, &test_config, lmp->version);

    delete lmp;
} 

TEST(FixLBMulticomponent, plain)
{
    LAMMPS *lmp = init_lammps();
    
    // Check for LAMMPS initailization with the fix, else skip the test
    if (!lmp) {                                                       
      std::cerr << "One or more prerequisite styles are not available "
                   "in this LAMMPS configuration:\n";               
      for (auto &prerequisite : test_config.prerequisites)         
        std::cerr << prerequisite.first << "_style " << prerequisite.second << "\n";
                                                                    
      GTEST_SKIP();                                                 
    }                                                                 
                                                                   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // get the base-class Fix * matching the style
    auto fixes = lmp->modify->get_fix_by_style("lb/multicomponent");    
        
    ASSERT_EQ(fixes.size(), 1);
    ASSERT_NE(fixes[0], nullptr);

    // cast to FixLbMulticomponent * is it exists
    auto *fix = dynamic_cast<FixLbMulticomponent *>(fixes[0]);
    
    // verify the cast to FixLbMulticomponent succeeded
    ASSERT_NE(fix, nullptr);
    
    // obtain initial total momentum
    double jx0, jy0, jz0, m0;
    jx0 = jy0 = jz0 = m0 = 0.0;
    fix->get_total_momentum(jx0, jy0, jz0);
    fix->get_total_mass(m0);
    if (rank == 0){
      std::cout << "Initial P: " << jx0 << " " << jy0 << " " << jz0 << std::endl;
      std::cout << "Initial M: " << m0 << std::endl;
      EXPECT_NEAR(jx0, 0, 1e-6);
      EXPECT_NEAR(jy0, 0, 1e-6);
      EXPECT_NEAR(jz0, 0, 1e-6);
    }
    
    // runs lammps time steps invoking fix_lb_multicomponent
    run_lammps(lmp); 

    // obtain final total momentum
    double jx1, jy1, jz1, m1;
    jx1 = jy1 = jz1 = m1 = 0.0;
    fix->get_total_momentum(jx1, jy1, jz1);
    fix->get_total_mass(m1);
    
    // verify conservation total momentum
    if (rank == 0){
      std::cout << "Final P: " << jx1 << " " << jy1 << " " << jz1 << std::endl;
      std::cout << "Final M: " << m1 << std::endl;
      EXPECT_NEAR(jx1, 0, 1e-6);
      EXPECT_NEAR(jy1, 0, 1e-6);
      EXPECT_NEAR(jz1, 0, 1e-6);
      EXPECT_NEAR(m0, m1, 1e-6);
    }
    delete lmp;
}

/*
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    MPI_Finalize();
    return result;
}
*/

