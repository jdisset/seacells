// Gaga: lightweight simple genetic algorithm library
// Copyright (c) Jean Disset 2015, All rights reserved.

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library.

#ifndef GAMULTI_HPP
#define GAMULTI_HPP

/****************************************
 *       TO ENABLE PARALLELISATION
 * *************************************/
// before including this file,
// #define OMP if you want OpenMP parallelisation
// #define CLUSTER if you want MPI parralelisation
#ifdef CLUSTER
#include <mpi.h>
#include <cstring>
#endif
#ifdef OMP
#include <omp.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <vector>
#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <random>
#include <utility>
#include <map>
#include <string>
#include "json/json.hpp"

#define PURPLE "\033[35m"
#define PURPLEBOLD "\033[1;35m"
#define BLUE "\033[34m"
#define BLUEBOLD "\033[1;34m"
#define GREY "\033[30m"
#define GREYBOLD "\033[1;30m"
#define YELLOW "\033[33m"
#define YELLOWBOLD "\033[1;33m"
#define RED "\033[31m"
#define REDBOLD "\033[1;31m"
#define CYAN "\033[36m"
#define CYANBOLD "\033[1;36m"
#define GREEN "\033[32m"
#define GREENBOLD "\033[1;32m"
#define NORMAL "\033[0m"

namespace GAGA {

using std::vector;
using std::string;
using std::unordered_set;
using std::map;
using std::unordered_map;
using std::cout;
using std::cerr;
using std::endl;
using fpType = std::vector<std::vector<double>>;  // footprints for novelty
using json = nlohmann::json;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;
using std::chrono::system_clock;

/******************************************************************************************
 *                                 GAGA LIBRARY
 *****************************************************************************************/
// This file contains :
// 1 - the Individual class template : an individual's generic representation, with its
// dna, fitnesses and
// behavior footprints (for novelty)

// 2 - the main GA class template
//
// About parallelisation :
// before including this file,
// #define OMP if you want OpenMP parallelisation
// #define CLUSTER if you want MPI parralelisation

/*****************************************************************************
 *                         INDIVIDUAL CLASS
 * **************************************************************************/
// A valid DNA class must have (see examples folder):
// DNA mutate()
// DNA crossover(DNA& other)
// static DNA random(int argc, char** argv)
// json& constructor
// void reset()
// json toJson()

template <typename DNA> struct Individual {
	DNA dna;
	map<string, double> fitnesses;  // map {"fitnessCriterName" -> "fitnessValue"}
	fpType footprint;               // individual's footprint for novelty computation
	string infos;                   // custom infos, description, whatever...
	bool evaluated = false;
	bool wasAlreadyEvaluated = false;
	double evalTime = 0.0;

	Individual() {}
	explicit Individual(const DNA &d) : dna(d) {}

	explicit Individual(const json &o) {
		assert(o.count("dna"));
		dna = DNA(o.at("dna").dump());
		if (o.count("footprint")) footprint = o.at("footprint").get<fpType>();
		if (o.count("fitnesses")) fitnesses = o.at("fitnesses").get<decltype(fitnesses)>();
		if (o.count("infos")) infos = o.at("infos");
		if (o.count("evaluated")) evaluated = o.at("evaluated");
		if (o.count("alreadyEval")) wasAlreadyEvaluated = o.at("alreadyEval");
		if (o.count("evalTime")) evalTime = o.at("evalTime");
	}

	// Exports individual to json
	json toJSON() const {
		json o;
		o["dna"] = json::parse(dna.toJSON());
		o["fitnesses"] = fitnesses;
		o["footprint"] = footprint;
		o["infos"] = infos;
		o["evaluated"] = evaluated;
		o["alreadyEval"] = wasAlreadyEvaluated;
		o["evalTime"] = evalTime;
		return o;
	}

	// Exports a vector of individual to json
	static json popToJSON(const vector<Individual<DNA>> &p) {
		json o;
		json popArray;
		for (auto &i : p) popArray.push_back(i.toJSON());
		o["population"] = popArray;
		return o;
	}

	// Loads a vector of individual from json
	static vector<Individual<DNA>> loadPopFromJSON(const json &o) {
		assert(o.count("population"));
		vector<Individual<DNA>> res;
		json popArray = o.at("population");
		for (auto &ind : popArray) res.push_back(Individual<DNA>(ind));
		return res;
	}
};

/*********************************************************************************
 *                                 GA CLASS
 ********************************************************************************/
// DNA requirements : see Individual class;
//
// Evaluaor class requirements (see examples folder):
// constructor(int argc, char** argv)
// void operator()(const Individual<DNA>& ind)
// const string name
//
// TYPICAL USAGE :
//
// GA<DNAType, EvalType> ga;
// ga.setPopSize(400);
// return ga.start();

template <typename DNA, typename Evaluator> class GA {
 protected:
	/*********************************************************************************
	 *                            MAIN GA SETTINGS
	 ********************************************************************************/
	bool novelty = false;  // is novelty enabled ?
	unsigned int verbosity =
	    2;  // 0 = silent; 1 = generations stats; 2 = individuals stats; 3 = everything
	unsigned int popSize = 500;       // nb of individuals in the population
	unsigned int nbElites = 1;        // nb of elites to keep accross generations
	unsigned int nbSavedElites = 1;   // nb of elites to save
	unsigned int tournamentSize = 3;  // nb of competitors in tournament
	unsigned int nbGen = 500;         // nb of generations
	double minNoveltyForArchive = 1;  // min novelty for being added to the general archive
	unsigned int KNN = 15;            // size of the neighbourhood for novelty
	bool savePopEnabled = true;       // save the whole population?
	bool saveArchiveEnabled = true;   // save the novelty archive?
	unsigned int saveInterval = 1;    // interval between 2 whole population saves
	string folder = "../evos/";       // where to save the results
	double crossoverProba = 0.2;      // crossover probability
	double mutationProba = 0.5;       // mutation probablility
	unordered_map<string, double>
	    proportions;  // {{"baseObj", 0.25}, {"novelty", 0.75}};  // fitness weight
	                  // proportions contains the relative weights of the objectives
	                  // if an objective is not present here but still used at evaluation
	                  // a default non weighted average will be used
	/********************************************************************************
	 *                                 SETTERS
	 ********************************************************************************/
 public:
	void enableNovelty() { novelty = true; }
	void disableNovelty() { novelty = false; }
	void enablePopulationSave() { savePopEnabled = true; }
	void disablePopulationSave() { savePopEnabled = false; }
	void enableArchiveSave() { saveArchiveEnabled = true; }
	void disableArchiveSave() { saveArchiveEnabled = false; }
	void setVerbosity(unsigned int lvl) { verbosity = lvl <= 3 ? (lvl >= 0 ? lvl : 0) : 3; }
	void setPopSize(unsigned int s) { popSize = s; }
	void setNbElites(unsigned int n) { nbElites = n; }
	void setNbSavedElites(unsigned int n) { nbSavedElites = n; }
	void setTournamentSize(unsigned int n) { tournamentSize = n; }
	void setNbGenerations(unsigned int n) { nbGen = n; }
	void setKNN(unsigned int n) { KNN = n; }
	void setPopSaveInterval(unsigned int n) { saveInterval = n; }
	void setSaveFolder(string s) { folder = s; }
	void setCrossoverProba(double p) {
		crossoverProba = p <= 1.0 ? (p >= 0.0 ? p : 0.0) : 1.0;
	}
	void setMutationProba(double p) {
		mutationProba = p <= 1.0 ? (p >= 0.0 ? p : 0.0) : 1.0;
	}
	void setMinNoveltyForArchive(double m) { minNoveltyForArchive = m; }
	void setObjectivesDistribution(map<string, double> d) { proportions = d; }
	void setObjectivesDistribution(string o, double d) { proportions[o] = d; }

	////////////////////////////////////////////////////////////////////////////////////

 protected:
	Evaluator evaluate;
	vector<Individual<DNA>>
	    archive;  // when novelty is enabled, we store the novel individuals there
	vector<Individual<DNA>> population;
	unsigned int currentGeneration = 0;
	// openmp/mpi stuff
	int procId = 0;
	int nbProcs = 1;
	int argc = 1;
	char **argv = nullptr;

	std::vector<std::map<std::string, std::map<std::string, double>>> genStats;

	std::random_device rd;
	std::default_random_engine globalRand = std::default_random_engine(rd());

	bool isBetter(double a, double b) { return a > b; }  // comparison btwn 2 fitnesses

 public:
	/*********************************************************************************
	 *                              CONSTRUCTOR
	 ********************************************************************************/
	GA(int ac, char **av) : evaluate(ac, av), argc(ac), argv(av) {
#ifdef CLUSTER
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
		MPI_Comm_rank(MPI_COMM_WORLD, &procId);
		if (procId == 0) {
			if (verbosity >= 3) {
				cout << "   -------------------" << endl;
				cout << CYAN << " MPI STARTED WITH " << NORMAL << nbProcs << CYAN << " PROCS "
				     << NORMAL << endl;
				cout << "   -------------------" << endl;
				cout << "Initialising population in master process" << endl;
			}
		}
#endif
	}

	/*********************************************************************************
	 *                          START THE BOUZIN
	 ********************************************************************************/
	// "Vroum vroum"
	int start() {
		bool finished = false;
		population.reserve(popSize);
		for (unsigned int i = 0; i < popSize; ++i) {
			population.push_back(Individual<DNA>(DNA::random(argc, argv)));
			population[population.size() - 1].evaluated = false;
		}
		if (procId == 0) {
			createFolder(folder);
			if (verbosity >= 1) printStart();
		}
		while (!finished) {
			auto tg0 = high_resolution_clock::now();
#ifdef CLUSTER
			MPI_distributePopulation();
#endif
#ifdef OMP
#pragma omp parallel for schedule(dynamic, 2)
#endif
			for (size_t i = 0; i < population.size(); ++i) {
				if (!population[i].evaluated) {
					auto t0 = high_resolution_clock::now();
					population[i].dna.reset();
					evaluate(population[i]);
					auto t1 = high_resolution_clock::now();
					population[i].evaluated = true;
					double indTime = std::chrono::duration<double>(t1 - t0).count();
					population[i].evalTime = indTime;
					population[i].wasAlreadyEvaluated = false;
				} else {
					population[i].evalTime = 0.0;
					population[i].wasAlreadyEvaluated = true;
				}
				if (verbosity >= 2) printIndividualStats(population[i]);
			}
#ifdef CLUSTER
			MPI_receivePopulation();
#endif
			if (procId == 0) {
				if (novelty) updateNovelty();
				auto tg1 = high_resolution_clock::now();
				double totalTime = std::chrono::duration<double>(tg1 - tg0).count();
				updateStats(totalTime);
				if (currentGeneration % saveInterval == 0) {
					if (savePopEnabled) savePop();
					if (novelty && saveArchiveEnabled) saveArchive();
				}
				if (verbosity >= 1) printGenStats(currentGeneration);
				saveBests(nbSavedElites);
				saveStats();
				prepareNextPop();
				finished = (currentGeneration++ >= nbGen);
			}
		}
#ifdef CLUSTER
		MPI_Finalize();
#endif
		return 0;
	}

// MPI specifics
#ifdef CLUSTER
	void MPI_distributePopulation() {
		if (procId == 0) {
			// if we're in the master process, we send b(i)atches to the others.
			// master will have the remaining
			unsigned int batchSize = population.size() / nbProcs;
			for (unsigned int dest = 1; dest < (unsigned int)nbProcs; ++dest) {
				vector<Individual<DNA>> batch;
				for (size_t ind = 0; ind < batchSize; ++ind) {
					batch.push_back(population.back());
					population.pop_back();
				}
				string batchStr = Individual<DNA>::popToJSON(batch).dump();
				std::vector<char> tmp(batchStr.begin(), batchStr.end());
				tmp.push_back('\0');
				MPI_Send(tmp.data(), tmp.size(), MPI_BYTE, dest, 0, MPI_COMM_WORLD);
			}
		} else {
			// we're in a slave process, we welcome our local population !
			int strLength;
			MPI_Status status;
			MPI_Probe(0, 0, MPI_COMM_WORLD, &status);  // we want to know its size
			MPI_Get_count(&status, MPI_CHAR, &strLength);
			char *popChar = new char[strLength + 1];
			MPI_Recv(popChar, strLength, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			// and we dejsonize !
			auto o = json::parse(popChar);
			population = Individual<DNA>::loadPopFromJSON(o);  // welcome bros!
			if (verbosity >= 3) {
				std::ostringstream buf;
				buf << endl
				    << "Proc " << PURPLE << procId << NORMAL << " : reception of "
				    << population.size() << " new individuals !" << endl;
				cout << buf.str();
			}
		}
	}

	void MPI_receivePopulation() {
		if (procId != 0) {  // if slave process we send our population to our mighty leader
			string batchStr = Individual<DNA>::popToJSON(population).dump();
			std::vector<char> tmp(batchStr.begin(), batchStr.end());
			tmp.push_back('\0');
			MPI_Send(tmp.data(), tmp.size(), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		} else {
			// master process receives all other batches
			for (unsigned int source = 1; source < (unsigned int)nbProcs; ++source) {
				int strLength;
				MPI_Status status;
				MPI_Probe(source, 0, MPI_COMM_WORLD, &status);  // determining batch size
				MPI_Get_count(&status, MPI_CHAR, &strLength);
				char *popChar = new char[strLength + 1];
				MPI_Recv(popChar, strLength + 1, MPI_BYTE, source, 0, MPI_COMM_WORLD,
				         MPI_STATUS_IGNORE);
				// and we dejsonize!
				auto o = json::parse(popChar);
				vector<Individual<DNA>> batch = Individual<DNA>::loadPopFromJSON(o);
				population.insert(population.end(), batch.begin(), batch.end());
				delete popChar;
				if (verbosity >= 3) {
					cout << endl
					     << "Proc " << procId << " : reception of " << batch.size()
					     << " treated individuals from proc " << source << endl;
				}
			}
		}
	}
#endif
	/*********************************************************************************
	 *                            NEXT POP GETTING READY
	 ********************************************************************************/
	// Là où qu'on fait les bébés.
	void prepareNextPop() {
		assert(population.size() > 0);
		assert(population.size() == popSize);
		vector<Individual<DNA>> nextGen;
		nextGen.reserve(popSize);
		vector<string> obj;
		for (auto &o : population[0].fitnesses) obj.push_back(o.first);
		auto subpop = multiObjTournament(obj, popSize);
		nextGen.insert(nextGen.end(), subpop.begin(), subpop.end());
		assert(nextGen.size() == popSize);
		population = nextGen;
	}

	void prepareSplitNextPop() {
		assert(population.size() > 0);
		assert(population.size() == popSize);
		// simplest MO
		vector<Individual<DNA>> nextGen;
		nextGen.reserve(popSize);
		vector<string> obj;
		for (auto &o : population[0].fitnesses) obj.push_back(o.first);
		size_t division = population.size() / obj.size();
		size_t remaining = population.size() % obj.size();
		for (size_t i = 0; i < obj.size(); ++i) {
			size_t extra = i < remaining ? 1 : 0;
			auto subpop = classicTournament(obj[i], division + extra);
			nextGen.insert(nextGen.end(), subpop.begin(), subpop.end());
		}
		assert(nextGen.size() == popSize);
		population = nextGen;
	}

	vector<Individual<DNA>> multiObjTournament(const std::vector<std::string> &objNames,
	                                           size_t n) {
		std::uniform_real_distribution<double> d(0.0, 1.0);
		assert(n > nbElites);
		vector<Individual<DNA>> newPop;
		newPop.reserve(n);
		// grab elites & add them to newPop
		auto elites = getElites(objNames, nbElites);
		for (auto &e : elites) {
			for (auto &i : e.second) {
				newPop.push_back(i);
			}
		}
		while (newPop.size() < n) {
			// choose 2 parents
			std::uniform_int_distribution<int> dint(0, population.size() - 1);
			std::vector<int> t0, t1;
			assert(tournamentSize > 0);
			for (unsigned int i = 0; i < tournamentSize; ++i) {
				t0.push_back(dint(globalRand));
				t1.push_back(dint(globalRand));
			}
			Individual<DNA> *p0 = &population[t0[0]];
			Individual<DNA> *p1 = &population[t1[0]];
			// we pick the objectives randomly
			std::uniform_int_distribution<int> dObj(0, objNames.size() - 1);
			std::string obj0 = objNames[dObj(globalRand)];
			std::string obj1 = objNames[dObj(globalRand)];
			for (unsigned int i = 0; i < tournamentSize; ++i) {
				if (isBetter(population[t0[i]].fitnesses.at(obj0), p0->fitnesses.at(obj0)))
					p0 = &population[t0[i]];
				if (isBetter(population[t1[i]].fitnesses.at(obj1), p1->fitnesses.at(obj1)))
					p1 = &population[t1[i]];
			}
			Individual<DNA> offspring;
			// create 1 offspring or simply copy one parent
			if (d(globalRand) < crossoverProba) {
				offspring = Individual<DNA>(p0->dna.crossover(p1->dna));
				offspring.evaluated = false;
			} else {
				offspring = *p0;
			}
			// mutate offspring
			if (d(globalRand) < mutationProba) {
				offspring.dna.mutate();
				offspring.evaluated = false;
			}
			newPop.push_back(offspring);
		}
		assert(newPop.size() == n);
		return newPop;
	}

	vector<Individual<DNA>> classicTournament(const std::string &objName, size_t n) {
		std::uniform_real_distribution<double> d(0.0, 1.0);
		assert(n > nbElites);
		vector<Individual<DNA>> newPop;
		newPop.reserve(n);
		// grab elites & add them to newPop
		auto elites = getElites({{objName}}, nbElites);
		for (auto &e : elites) {
			for (auto &i : e.second) {
				newPop.push_back(i);
			}
		}
		while (newPop.size() < n) {
			// choose 2 parents
			std::uniform_int_distribution<int> dint(0, population.size() - 1);
			std::vector<int> t0, t1;
			assert(tournamentSize > 0);
			for (unsigned int i = 0; i < tournamentSize; ++i) {
				t0.push_back(dint(globalRand));
				t1.push_back(dint(globalRand));
			}
			Individual<DNA> *p0 = &population[t0[0]];
			Individual<DNA> *p1 = &population[t1[0]];
			for (unsigned int i = 0; i < tournamentSize; ++i) {
				if (isBetter(population[t0[i]].fitnesses.at(objName), p0->fitnesses.at(objName)))
					p0 = &population[t0[i]];
				if (isBetter(population[t1[i]].fitnesses.at(objName), p1->fitnesses.at(objName)))
					p1 = &population[t1[i]];
			}
			Individual<DNA> offspring;
			// create 1 offspring or simply copy one parent
			if (d(globalRand) < crossoverProba) {
				offspring = Individual<DNA>(p0->dna.crossover(p1->dna));
				offspring.evaluated = false;
			} else {
				offspring = *p0;
			}
			// mutate offspring
			if (d(globalRand) < mutationProba) {
				offspring.dna.mutate();
				offspring.evaluated = false;
			}
			newPop.push_back(offspring);
		}
		assert(newPop.size() == n);
		return newPop;
	}

	unordered_map<string, vector<Individual<DNA>>> getElites(const vector<string> &obj,
	                                                         size_t n) {
		if (verbosity >= 3) {
			cerr << "getElites : nbObj = " << obj.size() << " n = " << n << endl;
		}
		unordered_map<string, vector<Individual<DNA>>> elites;
		for (auto &o : obj) {
			elites[o] = vector<Individual<DNA>>();
			elites[o].push_back(population[0]);
			size_t worst = 0;
			for (size_t i = 1; i < n && i < population.size(); ++i) {
				elites[o].push_back(population[i]);
				if (isBetter(elites[o][worst].fitnesses.at(o), population[i].fitnesses.at(o)))
					worst = i;
			}
			for (size_t i = n; i < population.size(); ++i) {
				if (isBetter(population[i].fitnesses.at(o), elites[o][worst].fitnesses.at(o))) {
					elites[o][worst] = population[i];
					for (size_t j = 0; j < n; ++j) {
						if (isBetter(elites[o][worst].fitnesses.at(o), elites[o][j].fitnesses.at(o)))
							worst = j;
					}
				}
			}
		}
		return elites;
	}

 protected:
	/*********************************************************************************
	 *                          NOVELTY RELATED METHODS
	 ********************************************************************************/
	// Novelty works with footprints. A footprint is just a vector of vector of doubles.
	// It is recommended that those doubles are within a same order of magnitude.
	// Each vector<double> is a "snapshot": it represents the state of the evaluation of
	// one individual at a certain time. Thus, a complete footprint is a combination
	// of one or more snapshot taken at different points in the
	// simulation (a vector<vector<double>>).
	// Snapshot must be of same size accross individuals.
	// Footprint must be set in the evaluator (see examples)

	static double getFootprintDistance(const fpType &f0, const fpType &f1) {
		assert(f0.size() == f1.size());
		double d = 0;
		for (size_t i = 0; i < f0.size(); ++i) {
			for (size_t j = 0; j < f0[i].size(); ++j) {
				d += std::pow(f0[i][j] - f1[i][j], 2);
			}
		}
		return sqrt(d);
	}

	// computeAvgDist (novelty related)
	// returns the average distance of a footprint fp to its k nearest neighbours
	// in an archive of footprints
	static double computeAvgDist(unsigned int K, const vector<Individual<DNA>> &arch,
	                             const fpType &fp) {
		double avgDist = 0;
		if (arch.size() > 1) {
			unsigned int k = arch.size() < K ? arch.size() : K;
			vector<Individual<DNA>> knn;
			knn.reserve(k);
			vector<double> knnDist;
			knnDist.reserve(k);
			std::pair<double, size_t> worstKnn = {getFootprintDistance(fp, arch[0].footprint),
			                                      0};  // maxKnn is the worst among the knn
			for (unsigned int i = 0; i < k; ++i) {
				knn.push_back(arch[i]);
				double d = getFootprintDistance(fp, arch[i].footprint);
				knnDist.push_back(d);
				if (d > worstKnn.first) {
					worstKnn = {d, i};
				}
			}
			for (size_t i = k; i < arch.size(); ++i) {
				double d = getFootprintDistance(fp, arch[i].footprint);
				if (d < worstKnn.first) {  // this one is closer than our worst knn
					knn[worstKnn.second] = arch[i];
					knnDist[worstKnn.second] = d;
					worstKnn.first = d;
					// we update maxKnn
					for (size_t j = 0; j < knn.size(); ++j) {
						if (knnDist[j] > worstKnn.first) {
							worstKnn = {knnDist[j], j};
						}
					}
				}
			}
			assert(knn.size() == k);
			for (size_t i = 0; i < knn.size(); ++i) {
				assert(getFootprintDistance(fp, knn[i].footprint) == knnDist[i]);
				avgDist += knnDist[i];
			}
			avgDist /= static_cast<double>(knn.size());
		}
		return avgDist;
	}
	void updateNovelty() {
		if (verbosity >= 2) {
			cout << endl
			     << endl;
			std::stringstream output;
			cout << GREY << " ❯❯  " << YELLOW << "COMPUTING NOVELTY " << NORMAL << " ⤵  "
			     << endl
			     << endl;
		}
		auto savedArchiveSize = archive.size();
		for (auto &ind : population) {
			archive.push_back(ind);
		}
		std::pair<Individual<DNA> *, double> best = {&population[0], 0};
		vector<Individual<DNA>> toBeAdded;
		for (auto &ind : population) {
			double avgD = computeAvgDist(KNN, archive, ind.footprint);
			bool added = false;
			if (avgD > minNoveltyForArchive) {
				toBeAdded.push_back(ind);
				added = true;
			}
			if (avgD > best.second) best = {&ind, avgD};
			if (verbosity >= 2) {
				std::stringstream output;
				output << GREY << " ❯ " << endl
				       << NORMAL << ind.infos << endl
				       << " -> Novelty = " << CYAN << avgD << GREY
				       << (added ? " (added to archive)" : " (too low for archive)") << NORMAL;
				if (verbosity >= 3)
					output << "Footprint was : " << footprintToString(ind.footprint);
				output << endl;
				std::cout << output.str();
			}
			ind.fitnesses["novelty"] = avgD;
		}
		archive.resize(savedArchiveSize);
		archive.insert(std::end(archive), std::begin(toBeAdded), std::end(toBeAdded));
		if (verbosity >= 2) {
			std::stringstream output;
			output << " Added " << toBeAdded.size() << " new footprints to the archive."
			       << std::endl
			       << "New archive size = " << archive.size() << " (was " << savedArchiveSize
			       << ")." << std::endl;
			std::cout << output.str() << std::endl;
		}
		if (verbosity >= 2) {
			std::stringstream output;
			output << "Most novel individual (novelty = " << best.second
			       << "): " << best.first->infos << endl;
			cout << output.str();
		}
	}

	// panpan cucul
	static inline string footprintToString(const vector<vector<double>> &f) {
		std::ostringstream res;
		res << "👣  " << json(f).dump();
		return res.str();
	}

	/*********************************************************************************
	 *                           STATS, LOGS & PRINTING
	 ********************************************************************************/
	void printStart() {
		int nbCol = 55;
		std::cout << std::endl
		          << GREY;
		for (int i = 0; i < nbCol - 1; ++i) std::cout << "━";
		std::cout << std::endl;
		std::cout << YELLOW << "              ☀     " << NORMAL << " Starting GAGA " << YELLOW
		          << "    ☀ " << NORMAL;
		std::cout << std::endl;
		std::cout << BLUE << "                      ¯\\_ಠ ᴥ ಠ_/¯" << std::endl
		          << GREY;
		for (int i = 0; i < nbCol - 1; ++i) std::cout << "┄";
		std::cout << std::endl
		          << NORMAL;
		std::cout << "  ▹ population size = " << BLUE << popSize << NORMAL << std::endl;
		std::cout << "  ▹ nb of elites = " << BLUE << nbElites << NORMAL << std::endl;
		std::cout << "  ▹ nb of tournament competitors = " << BLUE << tournamentSize << NORMAL
		          << std::endl;
		std::cout << "  ▹ mutation rate = " << BLUE << mutationProba << NORMAL << std::endl;
		std::cout << "  ▹ crossover rate = " << BLUE << crossoverProba << NORMAL << std::endl;
		std::cout << "  ▹ writing results in " << BLUE << folder << NORMAL << std::endl;
		if (novelty) {
			std::cout << "  ▹ novelty is " << GREEN << "enabled" << NORMAL << std::endl;
			std::cout << "    - KNN size = " << BLUE << KNN << NORMAL << std::endl;
		} else {
			std::cout << "  ▹ novelty is " << RED << "disabled" << NORMAL << std::endl;
		}
#ifdef CLUSTER
		std::cout << "  ▹ MPI parralelisation is " << GREEN << "enabled" << NORMAL
		          << std::endl;
#else
		std::cout << "  ▹ MPI parralelisation is " << RED << "disabled" << NORMAL
		          << std::endl;
#endif
#ifdef OMP
		std::cout << "  ▹ OpenMP parralelisation is " << GREEN << "enabled" << NORMAL
		          << std::endl;
#else
		std::cout << "  ▹ OpenMP parralelisation is " << RED << "disabled" << NORMAL
		          << std::endl;
#endif
		std::cout << GREY;
		for (int i = 0; i < nbCol - 1; ++i) std::cout << "━";
		std::cout << std::endl
		          << NORMAL;
	}
	void updateStats(double totalTime) {
		// stats organisations :
		// "global" -> {"genTotalTime", "indTotalTime", "maxTime", "nEvals", "nObjs"}
		// "obj_i" -> {"avg", "worst", "best"}
		assert(population.size());
		std::map<std::string, std::map<std::string, double>> currentGenStats;
		currentGenStats["global"]["genTotalTime"] = totalTime;
		double indTotalTime = 0.0, maxTime = 0.0;
		int nEvals = 0;
		int nObjs = population[0].fitnesses.size();
		for (const auto &o : population[0].fitnesses) {
			currentGenStats[o.first] = {
			    {{"avg", 0.0}, {"worst", o.second}, {"best", o.second}}};
		}
		for (const auto &ind : population) {
			indTotalTime += ind.evalTime;
			for (const auto &o : ind.fitnesses) {
				currentGenStats[o.first].at("avg") +=
				    (o.second / static_cast<double>(population.size()));
				if (isBetter(o.second, currentGenStats[o.first].at("best")))
					currentGenStats[o.first].at("best") = o.second;
				if (!isBetter(o.second, currentGenStats[o.first].at("worst")))
					currentGenStats[o.first].at("worst") = o.second;
			}
			if (ind.evalTime > maxTime) maxTime = ind.evalTime;
			if (!ind.wasAlreadyEvaluated) ++nEvals;
		}
		currentGenStats["global"]["indTotalTime"] = indTotalTime;
		currentGenStats["global"]["maxTime"] = maxTime;
		currentGenStats["global"]["nEvals"] = nEvals;
		currentGenStats["global"]["nObjs"] = nObjs;
		genStats.push_back(currentGenStats);
	}

	void printGenStats(unsigned int n) {
		const unsigned int l = 80;
		std::cout << tableHeader(l);
		std::ostringstream output;
		const auto &globalStats = genStats[n].at("global");
		output << "Generation " << CYANBOLD << n << NORMAL << " ended in " << BLUE
		       << globalStats.at("genTotalTime") << NORMAL << "s";
		std::cout << tableCenteredText(l, output.str(), BLUEBOLD NORMAL BLUE NORMAL);
		output = std::ostringstream();
		output << GREYBOLD << "(" << globalStats.at("nEvals") << " evaluations, "
		       << globalStats.at("nObjs") << " objs)" << NORMAL;
		std::cout << tableCenteredText(l, output.str(), GREYBOLD NORMAL);
		std::cout << tableSeparation(l);
		double timeRatio = 0;
		if (globalStats.at("genTotalTime") > 0)
			timeRatio = globalStats.at("indTotalTime") / globalStats.at("genTotalTime");
		output = std::ostringstream();
		output << "🕝  max: " << BLUE << globalStats.at("maxTime") << NORMAL << "s";
		output << ", 🕝  sum: " << BLUEBOLD << globalStats.at("indTotalTime") << NORMAL
		       << "s (x" << timeRatio << " ratio)";
		std::cout << tableCenteredText(l, output.str(), CYANBOLD NORMAL BLUE NORMAL "      ");
		std::cout << tableSeparation(l);
		for (const auto &o : genStats[n]) {
			if (o.first != "global") {
				output = std::ostringstream();
				output << GREYBOLD << "--◇" << GREENBOLD << std::setw(10) << o.first << GREYBOLD
				       << " ❯ " << NORMAL << " worst: " << YELLOW << std::setw(12)
				       << o.second.at("worst") << NORMAL << ", avg: " << YELLOWBOLD
				       << std::setw(12) << o.second.at("avg") << NORMAL << ", best: " << REDBOLD
				       << std::setw(12) << o.second.at("best") << NORMAL;
				std::cout << tableText(l, output.str(),
				                       "    " GREYBOLD GREENBOLD GREYBOLD NORMAL YELLOWBOLD NORMAL
				                           YELLOW NORMAL GREENBOLD NORMAL);
			}
		}
		std::cout << tableFooter(l);
	}

	void printIndividualStats(const Individual<DNA> &ind) {
		std::ostringstream output;
		output << GREYBOLD << "[" << YELLOW << procId << GREYBOLD << "]-▶ " << NORMAL;
		for (const auto &o : ind.fitnesses)
			output << " " << o.first << ": " << BLUEBOLD << std::setw(12) << o.second << NORMAL
			       << GREYBOLD << " |" << NORMAL;
		output << " 🕝 : " << BLUE << ind.evalTime << "s" << NORMAL;
		if (ind.wasAlreadyEvaluated)
			output << GREYBOLD << " (already evaluated)\n" << NORMAL;
		else
			output << "\n";
		if ((!novelty && verbosity >= 2) || verbosity >= 3) output << ind.infos << std::endl;
		std::cout << output.str();
	}

	std::string tableHeader(unsigned int l) {
		std::ostringstream output;
		output << "┌";
		for (auto i = 0u; i < l; ++i) output << "─";
		output << "┓\n";
		return output.str();
	}

	std::string tableFooter(unsigned int l) {
		std::ostringstream output;
		output << "┗";
		for (auto i = 0u; i < l; ++i) output << "─";
		output << "┛\n";
		return output.str();
	}

	std::string tableSeparation(unsigned int l) {
		std::ostringstream output;
		output << "|" << GREYBOLD;
		for (auto i = 0u; i < l; ++i) output << "-";
		output << NORMAL << "|\n";
		return output.str();
	}

	std::string tableText(int l, std::string txt, std::string unprinted = "") {
		std::ostringstream output;
		int txtsize = txt.size() - unprinted.size();
		output << "|" << txt;
		for (auto i = txtsize; i < l; ++i) output << " ";
		output << "|\n";
		return output.str();
	}

	std::string tableCenteredText(int l, std::string txt, std::string unprinted = "") {
		std::ostringstream output;
		int txtsize = txt.size() - unprinted.size();
		int space = l - txtsize;
		output << "|";
		for (int i = 0; i < space / 2; ++i) output << " ";
		output << txt;
		for (int i = (space / 2) + txtsize; i < l; ++i) output << " ";
		output << "|\n";
		return output.str();
	}

	/*********************************************************************************
	 *                         SAVING STUFF
	 ********************************************************************************/
	void saveBests(int n) {
		// save n bests dnas for all objectives
		vector<string> objectives;
		for (auto &o : population[0].fitnesses) {
			objectives.push_back(o.first);  // we need to know objective functions
		}
		auto elites = getElites(objectives, n);
		std::stringstream baseName;
		baseName << folder << "/gen" << currentGeneration;
		mkdir(baseName.str().c_str(), 0777);
		if (verbosity >= 3) {
			cerr << "created directory " << baseName.str() << endl;
		}
		for (auto &e : elites) {
			int id = 0;
			for (auto &i : e.second) {
				std::stringstream fileName;
				fileName << baseName.str() << "/" << e.first << "_" << i.fitnesses.at(e.first)
				         << "_" << id++ << ".dna";
				std::ofstream fs(fileName.str());
				if (!fs) {
					cerr << "Cannot open the output file." << endl;
				}
				fs << i.dna.toJSON();
				fs.close();
			}
		}
	}

	void saveStats() {
		std::stringstream csv;
		std::stringstream fileName;
		fileName << folder << "/stats.csv";
		csv << "generation";
		if (genStats.size() > 0) {
			for (auto &cat : genStats[0]) {
				std::stringstream column;
				column << cat.first << "_";
				for (auto &s : cat.second) {
					csv << "," << column.str() << s.first;
				}
			}
			csv << endl;
			for (size_t i = 0; i < genStats.size(); ++i) {
				csv << i;
				for (auto &cat : genStats[i]) {
					for (auto &s : cat.second) {
						csv << "," << s.second;
					}
				}
				csv << endl;
			}
		}
		std::ofstream fs(fileName.str());
		if (!fs) {
			cerr << "Cannot open the output file." << endl;
		}
		fs << csv.str();
		fs.close();
	}

	void createFolder(string baseFolder) {
		struct stat sb;
		mkdir(baseFolder.c_str(), 0777);
		auto now = system_clock::now();
		time_t now_c = system_clock::to_time_t(now);
		struct tm *parts = localtime(&now_c);

		std::stringstream fname;
		fname << evaluate.name << parts->tm_mday << "_" << parts->tm_mon + 1 << "_";
		int cpt = 0;
		std::stringstream ftot;
		do {
			ftot.clear();
			ftot.str("");
			ftot << baseFolder << fname.str() << cpt;
			cpt++;
		} while (stat(ftot.str().c_str(), &sb) == 0 && S_ISDIR(sb.st_mode));
		folder = ftot.str();
		mkdir(folder.c_str(), 0777);
	}

 public:
	void loadPop(string file) {
		std::ifstream t(file);
		std::stringstream buffer;
		buffer << t.rdbuf();
		auto o = json::parse(buffer.str());
		assert(o.count("population"));
		if (o.count("generation")) {
			currentGeneration = o.at("generation");
		} else {
			currentGeneration = 0;
		}
		population.clear();
		for (auto ind : o.at("population")) {
			population.push_back(Individual<DNA>(DNA(ind.at("dna"))));
			population[population.size() - 1].evaluated = false;
		}
	}

	void savePop() {
		json o = Individual<DNA>::popToJSON(population);
		o["evaluator"] = evaluate.name;
		o["generation"] = currentGeneration;
		std::stringstream baseName;
		baseName << folder << "/gen" << currentGeneration;
		mkdir(baseName.str().c_str(), 0777);
		std::stringstream fileName;
		fileName << baseName.str() << "/pop" << currentGeneration << ".pop";
		std::ofstream file;
		file.open(fileName.str());
		file << o.dump();
		file.close();
	}
	void saveArchive() {
		json o = Individual<DNA>::popToJSON(archive);
		o["evaluator"] = evaluate.name;
		std::stringstream baseName;
		baseName << folder << "/gen" << currentGeneration;
		mkdir(baseName.str().c_str(), 0777);
		std::stringstream fileName;
		fileName << baseName.str() << "/archive" << currentGeneration << ".pop";
		std::ofstream file;
		file.open(fileName.str());
		file << o.dump();
		file.close();
	}
};
}  // namespace GAGA
#endif
