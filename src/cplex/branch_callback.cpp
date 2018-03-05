/******************************************************************************
 * branch_callback.cpp: Implementation for Cplex Branch Callback.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jul 09, 2011 by andrade
 * Last update: Jul 13, 2011 by andrade
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

#include "branch_callback.hpp"

//------------------[ Default Constructor and Destructor ]-------------------//

BranchCallback::BranchCallback(IloEnv env, IloNumVarArray *ivars, std::vector<double> *ivals):
    IloCplex::BranchCallbackI(env), vars(ivars), vals(ivals), is_called(false)
    {}

BranchCallback::~BranchCallback() {}

//-------------------------[ Duplicate Callback ]----------------------------//

IloCplex::CallbackI* BranchCallback::duplicateCallback() const {
   return new (getEnv()) BranchCallback(*this);
}

//-----------------------[ main: called by Cplex ]---------------------------//

void BranchCallback::main() {
    #ifdef DEBUG
     cout << "\n--| BranchCallback" << endl;
    #endif

    vals->clear();

    for (int i = 0; i < vars->getSize(); ++i)
       vals->push_back(getValue((*vars)[i]));

    is_called = true;
    IloCplex::BranchCallbackI::abort();
}

//-------------------------------[ Reset ]-----------------------------------//

void BranchCallback::reset() {
    is_called = false;
}

