/******************************************************************************
 * take_duals.cpp: solve the linear model and take the dual variables for
 * each good.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 *  Created on : Mar 23, 2013 by andrade
 *  Last update: Mar 23, 2013 by andrade
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

#ifdef DEBUG
#include <iostream>
#include <algorithm>
#include <iterator>
#endif

#include <stdexcept>
#include <vector>
#include <cstring>

#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnostic pop

using namespace std;

//------------------------------[ Take duals ]--------------------------------//

void takeDualsValues(const char* filename, vector< double > &good_duals) {
    #ifdef DEBUG
    cout << "\n------------------------------------------------------\n"
         << "> Solving the linear programming " << filename << "\n"
         << "> and take the dual values."
         << endl;
    #endif

    IloEnv env;
    try {
        IloModel model(env);            // Problem model
        IloCplex cplex(env);            // Cplex algorithm
        IloNumVarArray variables(env);  // Variables from model
        IloObjective obj(env);          // Objective function
        IloRangeArray constraints(env); // Constraints

        #ifndef DEBUG
        // Suppress output
        cplex.setOut(env.getNullStream());
        cplex.setWarning(env.getNullStream());
        #else
        cplex.setParam(IloCplex::MIPDisplay, 4);
        #endif

        // Import/Build the model from LP file or CombinatorialAuction object.
        cplex.importModel(model, filename, obj, variables, constraints);

        // Relax all variables.
        model.add(IloConversion(env, variables, IloNumVar::Float));

        // Extract to algorithm.
        cplex.extract(model);

        // Turn off preprocessing.
        cplex.setParam(IloCplex::PreInd, false);
        cplex.setParam(IloCplex::Threads, 4);

        // Solve the relaxation.
        cplex.solve();

        IloNumArray duals(env);
        cplex.getDuals(duals, constraints);

        good_duals.clear();
        good_duals.reserve(duals.getSize());

        for(IloInt i = 0; i < duals.getSize(); ++i) {
            good_duals.push_back(duals[i]);
        }
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

    env.end();

    #ifdef DEBUG
    cout << "\n> Dual values:\n";
    copy(good_duals.begin(), good_duals.end(), ostream_iterator<double>(cout, " "));

    cout << "\n------------------------------------------------------"
         << endl;
    #endif
}
