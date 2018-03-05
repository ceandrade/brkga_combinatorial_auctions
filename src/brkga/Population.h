/*
 * Population.h
 *
 * Encapsulates a population of chromosomes represented by a vector of doubles. We don't decode
 * nor deal with random numbers here; instead, we provide private support methods to set the
 * fitness of a specific chromosome as well as access methods to each allele. Note that the BRKGA
 * class must have access to such methods and thus is a friend.
 *
 *  Created on : Jun 21, 2010 by rtoso
 *  Last update: Jun 22, 2011 by andrade
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 *      Collaborator: Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018
 * Rodrigo Franco Toso (rfrancotoso@gmail.com) and
 * Mauricio G.C. Resende
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef POPULATION_H
#define POPULATION_H

#include <vector>
#include <algorithm>
#include <exception>
#include <stdexcept>


class Population {
	template< class Decoder, class RNG, template<class T> class Compare >
	friend class BRKGA;

public:
	Population(const Population& other);
	Population(unsigned n, unsigned p);
	~Population();

	unsigned getN() const;	// Size of each chromosome
	unsigned getP() const;	// Size of population

	double operator()(unsigned i, unsigned j) const;	// Direct access to allele j of chromosome i

	// These methods REQUIRE fitness to be sorted, and thus a call to sortFitness() beforehand
	// (and I won't implement a lazy scheme to automatically update sortFitness() here).
	double getBestFitness() const;	// Returns the best fitness in this population
	double getFitness(unsigned i) const;	// Returns the fitness of chromosome i
	const std::vector< double >& getChromosome(unsigned i) const;	// Returns i-th best chromosome

private:
	std::vector< std::vector< double > > population;		// Population as vectors of prob.
	std::vector< std::pair< double, unsigned > > fitness;	// Fitness (double) of a each chromosome

	template< template<class T> class Compare >
	void sortFitness(); // Sorts 'fitness' by its first parameter by predicate template class.
	void setFitness(unsigned i, double f);				// Sets the fitness of chromosome i
	std::vector< double >& getChromosome(unsigned i);	// Returns a chromosome

	double& operator()(unsigned i, unsigned j);		// Direct access to allele j of chromosome i
	std::vector< double >& operator()(unsigned i);	// Direct access to chromosome i
};

Population::Population(const Population& pop) :
		population(pop.population),
		fitness(pop.fitness) {
}

Population::Population(const unsigned n, const unsigned p) :
		population(p, std::vector< double >(n, 0.0)), fitness(p) {
	if(p == 0) { throw std::range_error("Population size p cannot be zero."); }
	if(n == 0) { throw std::range_error("Chromosome size n cannot be zero."); }
}

Population::~Population() {
}

unsigned Population::getN() const {
	return population[0].size();
}

unsigned Population::getP() const {
	return population.size();
}

double Population::getBestFitness() const {
	return getFitness(0);
}

double Population::getFitness(unsigned i) const {
	return fitness[i].first;
}

const std::vector< double >& Population::getChromosome(unsigned i) const {
	return population[ fitness[i].second ];
}

std::vector< double >& Population::getChromosome(unsigned i) {
	return population[ fitness[i].second ];
}

void Population::setFitness(unsigned i, double f) {
	fitness[i].first = f;
	fitness[i].second = i;
}

template< template<class T> class Compare >
void Population::sortFitness() {
	sort(fitness.begin(), fitness.end(), Compare< std::pair< double, unsigned > >());
}

double Population::operator()(unsigned chromosome, unsigned allele) const {
	return population[chromosome][allele];
}

double& Population::operator()(unsigned chromosome, unsigned allele) {
	return population[chromosome][allele];
}

std::vector< double >& Population::operator()(unsigned chromosome) {
	return population[chromosome];
}

#endif
