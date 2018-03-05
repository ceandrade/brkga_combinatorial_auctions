/******************************************************************************
 * chromosome_generator.hpp: Interface for chromosome generator based on
 * LP relaxations of ILOG Cplex.
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

#ifndef CHROMOSOME_GENERATOR_HPP
#define CHROMOSOME_GENERATOR_HPP

#include <vector>

#pragma GCC diagnostic ignored "-Wignored-attributes"
#include <ilcplex/ilocplex.h>
#pragma GCC diagnostic pop

#include "combinatorial_auction.hpp"

/**
 * \brief Chromosome Generator based on LP relaxations.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2011
 *
 * This class generates chromosomes for BRKGA algorithm based on LP relaxations
 * of Ilog Cplex. It process only the first node of optimisation, putting cuts
 * as much as possible, and returning all solutions found (relaxated and
 * integer ones).
 */
class ChromosomeGenerator {
    public:
        /** \name Constructors and Destructor*/
        //@{
        /// Default constructor.
        ChromosomeGenerator();

        /// Destructor
        virtual ~ChromosomeGenerator();
        //@}

        /** \name Member Methods */
        //@{
        /** \brief Generate chromosomes from LP relaxations.
         * \param num_individuals max number of individuals to be generated.
         * \param time_limit overall generation time limit.
         * \param cplex_time_limit to Cplex try to find solutions.
         * \param cplex_iter_limit to Cplex try to generate cuts at root node.
         * \param pricing type of pricing used to generate individuals (first or second price)
         * \param auction to build the LP from it, in second price case.
         * \param filename The file that contains the problem description.
         * \return a reference to the generated chomosomes.
         */
        const std::vector< std::vector<double> >&
        generateChromosomes(unsigned num_individuals,
                            double time_limit,
                            double cplex_time_limit,
                            unsigned cplex_iter_limit,
                            Auction::Pricing pricing,
                            const CombinatorialAuction* auction = NULL,
                            const char* filename = NULL);
        //@}

    protected:
        /** \brief Build a second price model from combinatorial auction data
         * \param env Cplex environment
         * \param model Cplex model
         * \param vars the built variables
         * \param auction Auction data.
         */
        void buildSeconPriceModel(IloEnv env,
                                  IloModel model,
                                  IloNumVarArray vars,
                                  const CombinatorialAuction* auction);

        /// The generated chromosomes
        std::vector< std::vector<double> > chromosomes;
};

#endif // CHROMOSOME_GENERATOR_HPP

