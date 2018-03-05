/******************************************************************************
 * main.cpp: Main Process.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jan 20, 2014 by andrade
 * Last update: Jan 20, 2014 by andrade
 *
 * This code is released under LICENSE.md.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <iterator>
#include <boost/lexical_cast.hpp>

using namespace std;

#include "BRKGA.h"
#include "MTRand.h"
#include "timer.hpp"
#include "combinatorial_auction.hpp"
#include "combinatorial_auction_decoder.hpp"
#include "chromosome_generator.hpp"

void takeDualsValues(const char* filename, vector< double > &good_duals);

//------------------------[ Some control constants ]-------------------------//

// Controls stop criteria
enum StopRule{GENERATIONS = 'G', TARGET = 'T', IMPROVEMENT = 'I'};

// Initial population criteria
enum InitApproach{RAMDOM = 'R', LPRELAXATION = 'L'};

//-------------------------------[ Main ]------------------------------------//

int main(int argc, char* argv[])
{
//    unsigned num_individuals, double time_limit,
//    double cplex_time_limit, unsigned cplex_iter_limit, Auction::Pricing pricing,

    // First read parameters from command line:
    if(argc < 22) {
        cout << "usage: " << argv[0]
             << " <config-file> <seed> <stop-rule> <stop-arg> <decode-approach>"
             << " <initial-population> <init-percentage> <init-max-time> <init-cplex-time> <init-cplex-iters>"
             << " <max-time> <pricing-scheme> <instance-file> [instance-lp-file]"
             << "\nwhere: "
             << "\n - <config-file>: parameters of BRKGA algorithm"
             << "\n - <seed>: seed for random generator"
             << "\n - <stop-rule> <stop-arg>: stop rule and its arguments where"
             << "\n\t+ Generations <number_generations>: the algorithm runs until <number_generations>"
             << "\n\t+ Target <value of expected target>: runs until obtains the target value"
             << "\n\t+ Iterations <max generations without improvement>: runs until the solutions don't improve by max generations"
             << "\n - <decode-approach>: the approach used to make solutions:"
             << "\n\t+ Lazy: sort the bids in non-increasing order of they genes and in the following order"
             << "\n\t+ Greedy: sort the bids in non-increasing order of utility and choose the items in this order"
             << "\n\t+ Pirkul: sort the bids in non-increasing order of pseudo-utility(dual multipliers) and choose the items in this order"
             << "\n - <initial-population> <init-percentage> <init-iterations>: indicates how the algorithm must start the initial population"
             << "\n\t+ Random 0 0 0 0: generate all individuals randomly (the 0 has no effect)"
             << "\n\t+ LP <init-percentage> <init-max-time> <init-cplex-time> <init-cplex-iters>: use population-size <init-percentage> LP "
             << "\n\t  relaxations as start individuals (if not generate all, complete with random ones). Each call to Cplex is limited"
             << "\n\t  to <init-cplex-time>, that can uses <init-cplex-iters> tries to tight the model. This procedure are limited to"
             << "\n\t  <init-max-time> seconds. The algorithm uses <init-iterations> to tight the formulation and to generate each relaxation."
             << "\n - <max-time>: max running time (in seconds). If 0, the algorithm stops on chosen stop rule."
             << "\n - <pricing-scheme>: Princing scheme to the auction. Possible values are:"
             << "\n\t+ First: to first price scheme"
             << "\n\t+ Second: to second price scheme"
             << "\n - <instance-file>: instance file"
             << "\n - <instance-lp-file>: instance LP file (here, any format readable by Cplex). Mandatory only if initial population to be LP relaxations"

             << "\n - <pe>: % of elite items in the population"
             << "\n - <pm>: % of mutants introduced at each generation"
             << "\n - <rhoe>: prob. of inheriting each allele from elite parent"
             << "\n - <K>: number of independent populations"
             << "\n - <J>: interval at which elite chromosomes are exchanged (0 means no exchange)"
             << "\n - <M>: number of elite chromosomes to obtain from each population"
             << "\n - <E>: interval at which the populations are reset (0 means no reset)"
             << "\n\n ALL PARAMETERS ARE MANDATORY\n"
             << endl;
        return 64;  // BSD usage error code
    }

    // Loading parameters from command line
    const string configFile(argv[1]);       // config file
    unsigned long seed;                     // RNG seed
    char stop_rule;                         // stopping rule
    double stop_arg;                        // argument to stopping rule
    CombinatorialAuctionDecoder::DecoderApproach decode_approach; // true to greedy, false to lazy
    char init_population;                   // how initial population must be generated, ...
    double init_percentage;                 // ... how many individuals and...
    double init_max_time;                   // ... for how much time...
    double init_cplex_time;                 // ... using this time for each iteration...
    unsigned init_cplex_iters;              // ... how many iterations to be used
    double max_time;                        // max running time
    Auction::Pricing pricing;               // princing scheme
    const char* instance_file = argv[13];   // instance file
    const char* lp_file = argv[14];   // instance file

    try {
        seed = boost::lexical_cast<unsigned long>(argv[2]);
        stop_rule = StopRule(toupper(argv[3][0]));
        stop_arg = boost::lexical_cast<double>(argv[4]);
        decode_approach = CombinatorialAuctionDecoder::DecoderApproach(toupper(argv[5][0]));
        init_population = InitApproach(toupper(argv[6][0]));
        init_percentage = boost::lexical_cast<double>(argv[7]);
        init_max_time = boost::lexical_cast<double>(argv[8]);
        init_cplex_time = boost::lexical_cast<double>(argv[9]);
        init_cplex_iters = boost::lexical_cast<unsigned>(argv[10]);
        max_time = boost::lexical_cast<double>(argv[11]);
        pricing = Auction::Pricing(toupper(argv[12][0]));

        if((stop_rule != 'G' && stop_rule != 'T' && stop_rule != 'I') ||
           (decode_approach != CombinatorialAuctionDecoder::LAZY &&
            decode_approach != CombinatorialAuctionDecoder::GREEDY &&
            decode_approach != CombinatorialAuctionDecoder::PIRKUL) ||
           (init_population != 'R' && init_population != 'L') ||
           (pricing != Auction::FirstPrice && pricing != Auction::SecondPrice))

            throw exception();
    }
    catch(exception& e) {
        cout << "*** Bad arguments. Verify them!" << endl;
        return 64; // BSD usage error code
    }

    // Loading algorithm parameters from config file (code from rtoso).
    ifstream fin(configFile.c_str(), std::ios::in);
    if(!fin) {
        cout << "Cannot open configuration file: " << configFile << endl;
        return 66; // BSD file not found error code
    }
    else { fin.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit); }

    // Parameters read from configuration file:
    double pe;          // % of elite items in the population
    double pm;          // % of mutants introduced at each generation
    double rhoe;        // prob. of inheriting each allele from elite parent
    unsigned K;         // number of independent populations
    unsigned MAX_POP;   // maximum size of each population
    unsigned MAX_THR;   // number of threads in parallel decoding
    unsigned J;         // interval at which elite chromosomes are exchanged (0 means no exchange)
    unsigned M;         // number of elite chromosomes to obtain from each population
    unsigned E;         // interval at which the populations are reset (0 means no reset)

    try {
        string line;
        fin >> pe;      getline(fin, line); // read the parameter and ignore the rest of the line
        fin >> pm;      getline(fin, line);
        fin >> rhoe;    getline(fin, line);
        fin >> K;       getline(fin, line);
        fin >> MAX_POP; getline(fin, line);
        fin >> MAX_THR; getline(fin, line);
        fin >> J;       getline(fin, line);
        fin >> M;       getline(fin, line);
        fin >> E;

        fin.close();
    }
    catch(ifstream::failure& e) {
        cout << "Failure when reading configfile: " << e.what() << endl;
        fin.close();
        return 65;  // BSD file read error code
    }

    //-----------------------------------------//
    // Tunning

    pe = boost::lexical_cast<double>(argv[15]);
    pm = boost::lexical_cast<double>(argv[16]);
    rhoe = boost::lexical_cast<double>(argv[17]);
    K = boost::lexical_cast<unsigned>(argv[18]);
    J = boost::lexical_cast<unsigned>(argv[19]);
    M = boost::lexical_cast<unsigned>(argv[20]);
    E = boost::lexical_cast<unsigned>(argv[21]);

    //-----------------------------------------//

    //#ifdef DEBUG
    cout << "------------------------------------------------------"
         << "\n> Parameters\n"
         << "\nConfig file: " << configFile
         << "\n   + % of Elite: " << pe
         << "\n   + % of Mutants: " << pm
         << "\n   + Prob. inheritance (rhoe): " << rhoe
         << "\n   + # of independent populations: " << K
         << "\n   + maximum size of each population: " << MAX_POP
         << "\n   + # of threads: " << MAX_THR
         << "\n   + interval of chromossome exchange: " << J
         << "\n   + # of elite chromossome of each population: " << M
         << "\n   + reset interval: " << E
         << "\nSeed: " << seed
         << "\nStop Rule: [" << stop_rule << "] "
         << (stop_rule == 'G' ? "Generations -> " :
             (stop_rule == 'T' ? "Target -> " : "Improvement -> "))
         << stop_arg
         << "\nDecode approach: " << (decode_approach == CombinatorialAuctionDecoder::LAZY? "Lazy" :
                                     (decode_approach == CombinatorialAuctionDecoder::GREEDY? "Greedy" : "Pirkul"))
         << "\nInitial Population approach: " << (init_population == 'R'? "Random" : "LP Relaxations ");

    if(init_population == 'L')
        cout << "\n   + % of individuals: " << init_percentage
             << "\n   + Iterations: " << init_cplex_iters
             << "\n   + Init. Max time: " << init_max_time
             << "\n   + Init. Cplex time: " << init_cplex_time
             << "\n   + Init. Cplex iterarions: " << init_cplex_iters;

    cout  << "\nMax Time: " << max_time
          << "\nPricing Scheme: " << (pricing == Auction::FirstPrice? "First Price" : "Second Price")
          << "\n------------------------------------------------------" << endl;
    //#endif

    // Begin of real fun :-)
    try {
        // Loading instance file
        CombinatorialAuction auction;
        if(!auction.loadProblem(instance_file))
            return 65;  // BSD file read error code

        if(decode_approach == CombinatorialAuctionDecoder::PIRKUL)
            takeDualsValues(lp_file, auction.good_duals) ;

        // Decoder to auction
        CombinatorialAuctionDecoder decoder(&auction, pricing, decode_approach);
        // Random generator
        MTRand rng(seed);

        unsigned population_size = 10 * auction.num_bids;
        if(population_size > MAX_POP)
            population_size = (auction.num_bids > MAX_POP)? auction.num_bids : MAX_POP;

        // The algorithm
        BRKGA< CombinatorialAuctionDecoder, MTRand, greater >
            algorithm(auction.num_bids, population_size,
                      pe, pm, rhoe, decoder, rng, K, MAX_THR);

        // Hold best results
        double best_fitness = numeric_limits< double >::min();
        vector< double > best_chromosome;
        unsigned last_update_iteration = 0;
        double last_update_time = 0.0;
        unsigned update_offset = 0;
        unsigned large_offset = 0;
        double init_population_time = 0.0;

        Timer timer;

        // Use the LP relaxations as initial population
        if(init_population == LPRELAXATION) {
            if(argc < 10) {
                cout << "\n***********************************************************"
                     << "\n**** LP file not supplied"
                     << "\n***********************************************************"
                     << endl;
                return 65;  // BSD file read error code
            }

            try {
                ChromosomeGenerator cg;

                double perc = population_size * init_percentage;

                unsigned init_pop_size  =
                        (unsigned)(perc > 1.0 ? perc : 1.0);

                algorithm.setInitialPopulation(
                        cg.generateChromosomes(init_pop_size, init_max_time,
                                               init_cplex_time, init_cplex_iters,
                                               pricing, &auction, lp_file));
            }
            catch(exception& e) {
                cout << "\n***********************************************************"
                     << "\n****  Exception Occured: " << e.what()
                     << "\n***********************************************************"
                     << endl;
                return 70; // BSD software internal error code
            }
        }

        // Initialize properly and take elapsed time to do it
        algorithm.initialize();
        init_population_time = timer.elapsed();

        unsigned iteration = 1;
        bool run = true;

        #ifdef DEBUG
        cout << "\n------------------------------------------------------"
             << "\n>>>>> Evolving <<<<<" << endl;
        #endif

        cout << setiosflags(ios::fixed) << setprecision(6);

        // Reset the timer and do the iterations
        timer.restart();
        while(run) {
            // Run one iteration
            algorithm.evolve();

            // Elite-exchange:
            if(K > 0 && J > 0 && iteration % J == 0) {
                algorithm.exchangeElite(M);

                if(stop_rule != TARGET)
                    cout << "Exchanged " << M << " solutions from each population at iteration "
                         << iteration << "; best so far: " << algorithm.getBestFitness() << endl;
            }

            // Take the best solution
            double fitness = algorithm.getBestFitness();

            // Print some information
            #ifdef DEBUG
            cout << "Iteration " << iteration
                 << " / Best Current: " << fitness
                 << " / Best Overall: " << best_fitness << endl;
            #endif

            if((fitness - best_fitness) > 0.000001) {
                last_update_time = timer.elapsed();

                if(stop_rule != TARGET)
                    cout << "Improvement Iteration " << iteration
                         << " / Current Time: " << last_update_time
                         << " / Best Current: " << fitness
                         << " / Best Overall: " << best_fitness
                         << endl;

                best_fitness = fitness;
                best_chromosome = algorithm.getBestChromosome();
                update_offset = iteration - last_update_iteration;
                last_update_iteration = iteration;

                if(large_offset < update_offset)
                    large_offset = update_offset;
            }

            // Time to reset?
            unsigned iterWithoutImprovement = iteration - last_update_iteration;
            if(E > 0 && iterWithoutImprovement > 0 && iterWithoutImprovement % E == 0) {
                algorithm.reset();
                if(stop_rule != TARGET)
                    cout << ">> Reset chromosomes with random keys after " << E
                         << " iterations without improvement." << endl;
            }

            // Time to stop?
            switch(stop_rule) {
                case GENERATIONS:
                    if(iteration == (unsigned)stop_arg)
                        run = false;
                    break;

                case IMPROVEMENT:
                    if(iterWithoutImprovement >= (unsigned)stop_arg)
                        run = false;
                    break;

                case TARGET:
                    if(abs(stop_arg - best_fitness) < 0.000001)
                        run = false;
                    break;
            }

            if(timer.elapsed() > max_time) {
                run = false;
            }

            ++iteration;
        } // end of while
        #ifdef DEBUG
        cout << "\n\n> End of evolving"
             << "\n------------------------------------------------------" << endl;
        #endif

        // Take elapsed time to optmization
        double elapsed_time = timer.elapsed();

        if(stop_rule != TARGET) {
            // Print results
            list< unsigned > solution = decoder.getSolutionFromChromosome(best_chromosome);
            double value_by_auction = auction.getSolutionRevenue(solution, pricing);

            cout << "\nBest solution > value: " << value_by_auction << " | (bids): ";
            copy(solution.begin(), solution.end(), ostream_iterator< unsigned >(cout, " "));

            cout << "\n\nRevenue & Iterations & Last Update Iteration & Last Update Time & "
                 << "Update Offset & Large Offset & Init. Time (s) & Total Elap. Time (s)\n"
                 << setiosflags(ios::fixed) << setprecision(6)
                 << best_fitness << " & "
                 << setprecision(2)
                 << --iteration << " & "            // only to correct the iteration number
                 << last_update_iteration << " & "
                 << last_update_time << " & "
                 << update_offset << " & "
                 << large_offset << " & "
                 << init_population_time << " & "
                 << elapsed_time;
        }
        else {
            cout << setiosflags(ios::fixed) << setprecision(2)
                 << best_fitness << " & "
                 << setprecision(2)
                 << --iteration << " & "            // only to correct the iteration number
                 << init_population_time << " & "
                 << elapsed_time;
        }

        cout << "\n\n"
             << setiosflags(ios::fixed) << setprecision(6)
             << (-best_fitness);
        cout.flush();
    }
    catch(exception& e) {
        cout << "\n***********************************************************"
             << "\n****  Exception Occured: " << e.what()
             << "\n***********************************************************"
             << endl;
        return 70; // BSD software internal error code
    }
    return 0;
}

