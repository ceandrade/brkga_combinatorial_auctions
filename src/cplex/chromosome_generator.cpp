/******************************************************************************
 * chromosome_generator.cpp: Implementation for chromosome generator based on
 * LP relaxations of Ilog Cplex.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jul 11, 2011 by andrade
 * Last update: Jul 28, 2011 by andrade
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
#include <string>
#include <stdexcept>
#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
using std::cout;
using std::endl;
using std::string;

#ifdef DEBUG
#include <algorithm>
#include <iterator>
#include <limits>
#include <typeinfo>
using namespace std;
#endif

#include "chromosome_generator.hpp"
#include "branch_callback.hpp"
#include "heuristic_callback.hpp"
using std::vector;

//------------------[ Default Constructor and Destructor ]-------------------//

ChromosomeGenerator::ChromosomeGenerator(): chromosomes() {}
ChromosomeGenerator::~ChromosomeGenerator() {}

//-------------------------[ Generate Chromosomes ]--------------------------//

const std::vector< std::vector<double> >&
ChromosomeGenerator::generateChromosomes(unsigned num_individuals, double time_limit,
         double cplex_time_limit, unsigned cplex_iter_limit, Auction::Pricing pricing,
         const CombinatorialAuction* auction, const char* filename) {

    #ifdef DEBUG
     std::cout << "\n------------------------------------------------------\n"
               << "> Generating Chromosomes using "
               << (pricing == Auction::FirstPrice? "First Price" : "Second Price")
               << "\n";
     std::cout.flush();
    #endif

    // Clear previous chromosomes
    chromosomes.clear();

    // Some checks
    if(pricing == Auction::FirstPrice && filename == NULL) {
        throw std::runtime_error("The LP file hasn't been supplied");
    }

    // Some checks
    if(pricing == Auction::SecondPrice && auction == NULL) {
        throw std::runtime_error("The auction object hasn't been supplied");
    }

    // Cplex environment
    IloEnv env;
    try {
        IloModel model(env);            // Problem model
        IloCplex cplex(env);            // Cplex algorithm
        IloNumVarArray variables(env);  // Variables from model
        vector<double> values;          // Variables' values

        #ifndef DEBUG
        // Suppress output
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        #endif

        // Import/Build the model from LP file or CombinatorialAuction object
        if(pricing == Auction::FirstPrice) {
            IloObjective obj(env);      // Objective function
            IloRangeArray rngs(env);    // Row ranges
            cplex.importModel(model, filename, obj, variables, rngs);
        }
        // Second price: build and take only x variables
        else {
            IloNumVarArray temp_vars(env);  // Variables from model

            buildSeconPriceModel(env, model, temp_vars, auction);

            // Take x variables
            for(int i = 0; i < temp_vars.getSize(); ++i) {
              if(temp_vars[i].getName()[0] == 'x')
                  variables.add(temp_vars[i]);
            }
        }

        //extract to algorithm
        cplex.extract(model);
        //cplex.exportModel("teste.lp");
        //
        values.reserve(variables.getSize());

        // Create and set the heuristic and branch callbacks
        BranchCallback* branch_cb = new (env) BranchCallback(env, &variables, &values);
        HeuristicCallback* heuristic_cb = new (env) HeuristicCallback(env, &variables, &values, cplex_iter_limit, cplex_time_limit);

        cplex.use(branch_cb);
        cplex.use(heuristic_cb);

        // Turn off probing and set the number of threads
        cplex.setParam(IloCplex::Probe, -1);
        cplex.setParam(IloCplex::Threads, 4);

        // Stop at root node and look for good bounds
        cplex.setParam(IloCplex::NodeLim, 1);

        // Timer
        boost::timer timer;

        // We try to solve the root problem, fixing variables.
        int bound = 0;
        unsigned j = 0;
        for(int i = -1; (i < variables.getSize()) &&
                        (j < num_individuals) &&
                        (timer.elapsed() < time_limit) ; ) {

            // Try to fix variables. When i = -1, nothing is done.
            // After, we fix each variable to 0.0 or 1.0, alternating.
            if(i >= 0) {
                variables[i].setBounds(bound, bound);

                if(bound == 0) {
                    bound = 1;

                    // Restore the bounds of last variable
                    if(i > 0)
                        variables[i-1].setBounds(0.0, 1.0);
                }
                else {
                    bound = 0;
                    ++i;
                }
            }
            else { ++i; }

            #ifdef DEBUG
            cout << "\n\n------------------| Generating a new |------------------";

            if(i > 0) {
                cout << "\n>>> Variable " << i-1 << " ("  << variables[i-1].getLB()
                     << ", " << variables[i-1].getUB() << ")"
                     << endl << endl;
            }
            else
                cout << "\n>>> Free variables" << endl << endl;
            #endif

            // Reset callbacks
            branch_cb->reset();
            heuristic_cb->reset();

            // Solve
            if(cplex.solve()) {
                // If we have a relaxation, take it.
                if(branch_cb->is_called || heuristic_cb->is_called) {
                    chromosomes.push_back(values);
                    ++j;
                }

                // Take the best integer solution too
                if(cplex.getStatus() == IloAlgorithm::Feasible ||
                   cplex.getStatus() == IloAlgorithm::Optimal) {

                    values.clear();

                    double check = 0.0;
                    for(int i = 0; i < variables.getSize(); ++i) {
                       double v = cplex.getValue(variables[i]);
                       check += v;
                       values.push_back(v);
                    }

                    if(check > 0.000001) {
                        chromosomes.push_back(values);
                        ++j;
                    }
                }
            }

            #ifdef DEBUG
            cout << "\n\n++++++++++++ Values empty: " << values.empty()
                 << "\n++++++++++++ Status: " << cplex.getStatus()
                 << "\n++++++++++++ Cplex Status: " << cplex.getCplexStatus()
                 << "\n++++++++++++ Cplex getObjValue: " << cplex.getObjValue()
                 << "\n++++++++++++ Cplex getBestObjValue: " << cplex.getBestObjValue()
                 << "\n\n*********** Generated: " << chromosomes.size()
                 << "\n*********** Elapsed time: " << timer.elapsed()
                 << endl << endl;
            #endif
        } // endfor

        // We mustn't forget to finalize the cplex environment
        env.end();
    }
    catch(IloException& e) {
       cout << "\n***********************************************************"
            << "\n**** "  << e
            << "***********************************************************"
            << endl;

       // We mustn't forget to finalize the cplex environment
       env.end();
       throw;
    }

    #ifdef DEBUG
    cout << "\n> Generated: " << chromosomes.size()
         << "\n------------------------------------------------------" << endl;
    #endif

    return chromosomes;
}

//-----------------------[ Generate by Second Price ]------------------------//

void ChromosomeGenerator::buildSeconPriceModel(IloEnv env, IloModel model, IloNumVarArray vars,
                                               const CombinatorialAuction* auction) {

     // Clean the model and vars
     model.end();
     vars.endElements();

     model = IloModel(env);

     // The epsilon second price
     const double epsilon = auction->worst_offer / auction->num_bids;

     // Variables that represents whether a bid wins or not.
     IloNumVarArray var_bids(env);
     for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
         IloBoolVar var = IloBoolVar(env, (string("x") + boost::lexical_cast<std::string>(bid_j)).c_str());
         var_bids.add(var);
         vars.add(var);
     }

     // Variables that indicates if a bid i, called by candidate, intercepts bid j
     // can be used to calculate the payments of bid j.
     // To build these variables, we take the edges from intersection graph.
     IloArray< IloNumVarArray > var_candidates(env);
     for(unsigned i = 0; i < auction->num_bids; ++i) {
         var_candidates.add(IloNumVarArray(env, auction->num_bids));
     }

     for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
         for(vector< unsigned >::const_iterator
             it_bid_i = auction->edges[bid_j].begin();
             it_bid_i != auction->edges[bid_j].end(); ++it_bid_i) {

             IloBoolVar var =
                 IloBoolVar(env, (string("c") +
                                  boost::lexical_cast<std::string>(*it_bid_i) + "_" +
                                  boost::lexical_cast<std::string>(bid_j)).c_str());

             var_candidates[bid_j][*it_bid_i] = var;
             vars.add(var);
         }
     }

     // Variables that represents the payments to be done by winner bidders.
     // Payments are greater than or equal to zero, to ensure
     // "Individual Rationality" propriety, and less than or equal to the second
     // price bid related to the actual bid, to try invoke the "incentive-compatibility"
     // propriety.
     IloNumVarArray var_payments(env);

     // This variables are created to force the model to choose bids that pays nothing.
     // We create a epsilon second price to these bids.
     std::map< unsigned, IloBoolVar > var_dummy_candidates;

     for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
         double ub = epsilon;

         // Take the second price...
         if(auction->next_lower_bids[bid_j] != auction->edges[bid_j].end()) {
             ub = auction->bids[*(auction->next_lower_bids[bid_j])].value;
         }
         //... or the epsilon price. In this case, we need a dummy candidate.
         else {
             IloBoolVar var = IloBoolVar(env, (string("c_dummy") +
                                               boost::lexical_cast<std::string>(bid_j)).c_str());

             var_dummy_candidates[bid_j] = var;
             vars.add(var);
         }

         IloNumVar var = IloNumVar(env, 0, ub, IloNumVar::Float,
                                   (string("p") + boost::lexical_cast<std::string>(bid_j)).c_str());
         var_payments.add(var);
         vars.add(var);
     }

     // Constraints 1: all payments must be less or equal the respective
     // bids, i.e., generalized second prices. In fact, the code above
     // try to do it.
 //    for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
 //        model.add(var_payments[bid_j] <= auction->bids[bid_j].value * var_bids[bid_j]);
 //    }

     // Constraints 2: Ensure the bid i, candidate to calculate bid j payment,
     // hasn't been used to calculate the payment for other bid k.
     //TODO: This constraints aren't very good.
 //    for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
 //        for(vector< unsigned >::const_iterator
 //            it_bid_i = auction->edges[bid_j].begin();
 //            it_bid_i != auction->edges[bid_j].end(); ++it_bid_i) {
 //
 //            for(vector< unsigned >::const_iterator
 //                it_bid_k = auction->edges[*it_bid_i].begin();
 //                it_bid_k != auction->edges[*it_bid_i].end(); ++it_bid_k) {
 //
 //                if(*it_bid_k != bid_j) {
 //                    model.add(var_candidates[bid_j][*it_bid_i] <= (1 - var_bids[*it_bid_k]));
 //                }
 //            }
 //        }
 //    }

     // Constraints 3: take the payment from candidate
     for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
         IloExpr sum_candidates(env);

         // Sum all candidates for bid j
         for(vector< unsigned >::const_iterator
             it_bid_i = auction->edges[bid_j].begin();
             it_bid_i != auction->edges[bid_j].end(); ++it_bid_i) {

             sum_candidates += auction->bids[*it_bid_i].value * var_candidates[bid_j][*it_bid_i];
         }

         // If we have a dummy candidate, put it in the constraint
         if(var_dummy_candidates.find(bid_j) != var_dummy_candidates.end())
             sum_candidates += epsilon * var_dummy_candidates[bid_j];

         model.add(var_payments[bid_j] == sum_candidates);
         sum_candidates.end();
     }

     // Constraints 4: only one candidate must be chosen for each bid.
     // Constraints 4': if a bid is chosen, we need a candidate too.
     for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
         IloExpr sum_candidates(env);

         for(vector< unsigned >::const_iterator
             it_bid_i = auction->edges[bid_j].begin();
             it_bid_i != auction->edges[bid_j].end(); ++it_bid_i) {

             sum_candidates += var_candidates[bid_j][*it_bid_i];
         }

         // If we have a dummy candidate, put it in the constraint
         if(var_dummy_candidates.find(bid_j) != var_dummy_candidates.end())
             sum_candidates += epsilon * var_dummy_candidates[bid_j];

         model.add(sum_candidates <= 1);
         model.add(var_bids[bid_j] >= sum_candidates);

         sum_candidates.end();
     }

     // Constraints 5: stable set constraints.
     for(unsigned good = 0; good < auction->num_goods + auction->num_dummy_goods; ++good) {
         IloExpr sum_bids(env);

         for(list< unsigned >::const_iterator it_bid = auction->goods_bids[good].begin();
             it_bid != auction->goods_bids[good].end(); ++it_bid) {
             sum_bids += var_bids[*it_bid];
         }

         model.add(sum_bids <= 1);

         sum_bids.end();
     }

     //TODO: verify this constraint
     // Constraint 6: force that, at least, one bid be chosen
     //model.add(IloSum(var_bids) >= 1);

     //TODO: we can use these edge-by-edge constraints too.
     // To prevent repeated edges.
     //    map< unsigned, map< unsigned, bool > > repeat;
 //    for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j) {
 //        for(vector< unsigned >::const_iterator
 //            it_bid_i = auction->edges[bid_j].begin();
 //            it_bid_i != auction->edges[bid_j].end(); ++it_bid_i) {
 //
 //            if(!repeat[bid_j][*it_bid_i]) {
 //                model.add(var_bids[bid_j] + var_bids[*it_bid_i] <= 1);
 //
 //                repeat[bid_j][*it_bid_i] = true;
 //                repeat[*it_bid_i][bid_j] = true;
 //            }
 //        }
 //     }

     IloNumArray bid_values(env);
     for(unsigned bid_j = 0; bid_j < auction->num_bids; ++bid_j)
         bid_values.add(auction->bids[bid_j].value);

     // Objective function: maximize the payments.
     model.add(IloMaximize(env, IloSum(var_payments)));
     //model.add(IloMaximize(env, IloSum(var_payments) + IloSum(var_bids)));
     //model.add(IloMaximize(env, IloScalProd(bid_values, var_bids) + IloSum(var_payments)));
     //model.add(IloMaximize(env, IloScalProd(bid_values, var_bids)));
 }

