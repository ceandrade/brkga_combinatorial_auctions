/******************************************************************************
 * heuristic_callback.cpp: Implementation for Cplex Heuristic Callback.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jul 13, 2011 by andrade
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

#ifdef DEBUG
#include <iostream>
#include <algorithm>
#include <iterator>
#include <limits>
#include <stdexcept>
using namespace std;
#endif

#include "heuristic_callback.hpp"

//------------------[ Default Constructor and Destructor ]-------------------//

HeuristicCallback::HeuristicCallback(IloEnv env, IloNumVarArray *ivars,
                                     std::vector<double> *ivals, unsigned _max_iterations,
                                     double _max_time):
    IloCplex::HeuristicCallbackI(env), vars(ivars), vals(ivals),
    max_iterations(_max_iterations), max_time(_max_time), is_called(false),
    iteration(0), timer()
    {}

HeuristicCallback::~HeuristicCallback() {}

//-------------------------[ Duplicate Callback ]----------------------------//

IloCplex::CallbackI* HeuristicCallback::duplicateCallback() const {
   return new (getEnv()) HeuristicCallback(*this);
}

//-----------------------[ main: called by cplex ]---------------------------//

void HeuristicCallback::main() {
    #ifdef DEBUG
     cout << "\n--| HeuristicCallback> iterarion: " <<  iteration
          << " | elapsed time: " << timer.elapsed()
          << " | max time: " << max_time
          << endl;
    #endif

    if((++iteration >= max_iterations) || (timer.elapsed() >= max_time)) {
        vals->clear();
        for (int i = 0; i < vars->getSize(); ++i)
            vals->push_back(getValue((*vars)[i]));

        is_called = true;

        IloCplex::HeuristicCallbackI::abort();
    }
}

//-------------------------------[ Reset ]-----------------------------------//

void HeuristicCallback::reset() {
    is_called = false;
    iteration = 0;
    timer.restart();
}
