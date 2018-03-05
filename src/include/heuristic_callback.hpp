/******************************************************************************
 * heuristic_callback.hpp: Interface for Cplex Heuristic Callback.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2011-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jul 13, 2011 by andrade
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

#ifndef HEURISTIC_CALLBACK_HPP
#define HEURISTIC_CALLBACK_HPP

#include <vector>
#include <cstring>
#include <boost/timer.hpp>

#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnostic pop

/**
 * \brief Heuristic callback for Cplex.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2011
 *
 * This class contains a heuristic callback for Ilog Cplex. This callback extracts
 * information about LP relaxation of a standard combinatorial auction model
 * and translates these results on vectors of doubles.
 */
class HeuristicCallback: public IloCplex::HeuristicCallbackI {
    public:
        /** \name Constructors and Destructor */
        //@{
        /** Default constructor.
         * \param env Cplex environment
         * \param ivars Reference to model variables
         * \param ivals Reference to a vector that will hold the variables' values
         * \param _max_iterations Number of iterations to run.
         * \param _max_time Max time to run
         */
        HeuristicCallback(IloEnv env, IloNumVarArray *ivars,
                          std::vector<double> *ivals,
                          unsigned _max_iterations,
                          double _max_time);

        /// Destructor
        virtual ~HeuristicCallback();
        //@}

        /** \name Mandatory Cplex Methods */
        //@{
        /// Duplicate the callback
        IloCplex::CallbackI* duplicateCallback() const;

        /// Main method called by Cplex. Read values of root relaxation, then stop.
        virtual void main();
        //@}

        /** \name Member Methods */
        //@{
        /// Reset the counters and flags
        void reset();
        //@}

    public:
        IloNumVarArray *const vars;         ///< Variables from model
        std::vector<double> *const vals;    ///< Variables' values
        unsigned max_iterations;            ///< Number of iterations to run
        double max_time;                    ///< Max time to run
        bool is_called;                     ///< Indicates if this callback was called

    protected:
        unsigned iteration;   ///< Current iterarion

        boost::timer timer;   ///< timer
};

#endif // HEURISTIC_CALLBACK_HPP

