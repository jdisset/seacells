#ifndef SCENARIO_HPP
#define SCENARIO_HPP

#include "../external/cxxopts.hpp"
#include "typesconfig.hpp"
#include "config.hpp"
#include <mecacell/mecacell.h>
#include <mecacell/grid.hpp>
#include <chrono>
#include <sstream>
#include <string>

template <typename Cell> class Scenario {
	struct PosIntegrator {
		static void updatePosition(Cell& c, const float_t& dt) {
			// c.setVelocity(c.getVelocity() + c.getForce() * dt / c.getMass());
			c.setVelocity(MecaCell::Vec::zero());
			c.setPrevposition(c.getPosition());
			c.setPosition(c.getPosition() + (c.getForce() * dt * dt / c.getMass()));
		}

		static void updateOrientation(Cell& c, const float_t& dt) {
			// c.setAngularVelocity(c.getAngularVelocity() +
			// c.getTorque() * dt / c.getMomentOfInertia());
			c.setAngularVelocity(MecaCell::Vec::zero());
			c.setOrientationRotation(c.getOrientationRotation() +
			                         (c.getTorque() * dt * dt / c.getMomentOfInertia()));
			c.updateCurrentOrientation();
		}
	};

	struct NutrientSource {
		MecaCell::Vec pos;
		double sqradius;
		double content;
		double initialcontent = Config::NUTRIENT_QUANTITY;
	};

 public:
	using World = MecaCell::BasicWorld<Cell>;
	using CellType = Cell;
	using CtrlType = typename Cell::CtrlType;

 protected:
	World w;
	double simDuration = Config::DEFAULT_SIM_DURATION;
	unsigned int maxCells = Config::DEFAULT_MAX_CELLS;
	unique_ptr<Cell> stemCell;
	std::chrono::time_point<std::chrono::system_clock> start;
	int randomSeed = 1000;
	MecaCell::Vec stemCellPosition{0, 30, 0};
	MecaCell::Grid<Cell*> cellgrid =
	    MecaCell::Grid<Cell*>(2.0 * MecaCell::DEFAULT_CELL_RADIUS);

 public:
	std::vector<NutrientSource> nutrientSources;
	double simTime = 0.0;
	double plantEnergy = 0.0;
	unsigned int getMaxUpdates() { return simDuration / w.getDt(); }
	void setStemCell(Cell* c) { stemCell = unique_ptr<Cell>(c); }

	void init(int argc, char** argv) {
		start = std::chrono::system_clock::now();
		try {
			cxxopts::Options options(argv[0]);
			options.add_options()("duration", "simulation duration",
			                      cxxopts::value<double>(simDuration));
			options.add_options()("maxcell", "max number of cells",
			                      cxxopts::value<unsigned int>(maxCells));
			options.add_options()("r,random", "random grn as stem cell");
			options.add_options()("f,file", "stem cell dna from file",
			                      cxxopts::value<std::string>());
			options.parse(argc, argv);
			if (options.count("random")) {
				stemCell = unique_ptr<Cell>(new Cell(CtrlType::random(argc, argv)));
			} else if (options.count("file")) {
				std::ifstream fstr(options["file"].as<std::string>());
				std::stringstream buffer;
				buffer << fstr.rdbuf();
				stemCell = unique_ptr<Cell>(new Cell(CtrlType(buffer.str())));
			}
		} catch (const cxxopts::OptionException& e) {
			std::cout << "error parsing options: " << e.what() << std::endl;
			exit(1);
		} catch (const std::bad_cast& e) {
			std::cout << "bad cast: " << e.what() << std::endl;
			exit(1);
		}
		w.setDt(Config::SIM_DT);
		if (stemCell) {
			w.setViscosityCoef(0.0);  // we handle viscosity by ourselves
			stemCell->setPosition(MecaCell::Vec(0, Config::STEMCELL_Y, 0));
			stemCell->nutrientLevel[WATER] = Config::STEMCELL_NUT0;
			stemCell->nutrientLevel[LIGHT] = Config::STEMCELL_NUT1;
			w.addCell(new Cell(*stemCell));
			initNutrientsSources();
		} else {
			std::cerr << "No stem cell, aborting." << std::endl;
			exit(1);
		}
	}

	void initNutrientsSources() {
		nutrientSources = std::vector<NutrientSource>();
		if (Config::NB_NUTRIENTS_SOURCES > 0) {
			auto internalRand = std::mt19937(randomSeed);
			auto uniformDist = std::uniform_real_distribution<double>(-0.5, 0.5);
			nutrientSources.push_back({stemCell->getPosition(),
			                           std::pow(Config::TYPICAL_NUTRIENTS_RADIUS, 2),
			                           Config::NUTRIENT_QUANTITY, Config::NUTRIENT_QUANTITY});

			for (auto i = 1u; i < Config::NB_NUTRIENTS_SOURCES; ++i) {
				MecaCell::Vec boundingCubeCenter(0.0, -Config::FIRST_NUTRIENT_SOURCE_THRESHOLD -
				                                          Config::NUTRIENTS_BOUNDING_DEPTH * 0.5,
				                                 0.0);
				MecaCell::Vec rdmPos =
				    boundingCubeCenter +
				    MecaCell::Vec(uniformDist(internalRand) * Config::NUTRIENTS_BOUNDING_AREA,
				                  uniformDist(internalRand) * Config::NUTRIENTS_BOUNDING_DEPTH,
				                  uniformDist(internalRand) * Config::NUTRIENTS_BOUNDING_AREA);
				double qtty =
				    Config::NUTRIENT_QUANTITY *
				    (1.0 + (pow(abs(rdmPos.y() + Config::FIRST_NUTRIENT_SOURCE_THRESHOLD),
				                Config::NUTRIENT_DEPTH_INCREASE_POW) *
				            Config::NUTRIENT_DEPTH_INCREASE_COEF));
				nutrientSources.push_back(
				    {rdmPos, std::pow(Config::TYPICAL_NUTRIENTS_RADIUS, 2), qtty, qtty});
			}
		}
		// for (auto& n : nutrientSources) {
		// stemCell->setPosition(n.pos);
		// w.addCell(new Cell(*stemCell));
		//}
	}

	void updateCellsSensedNutrients() {
		for (auto& c : w.cells)
			c->setSensedNutrients(WATER, computeNutrientIntensity(c->getPosition()));
	}

	void terminate() {
		// auto end = std::chrono::system_clock::now();
		// std::chrono::duration<double> diff = end - start;
	}

	inline double diffusedQtty(double deltaP, double area, double l, double dt) {
		return -dt * (static_cast<double>(Config::NUTRIENTS_DIFFUSION_K) * area * deltaP) /
		       std::max(1.0, (static_cast<double>(Config::NUTRIENTS_VISCOSITY) *
		                      std::max(1.0, l)));
	}

	void shineOn() {
		const double resolution = MecaCell::DEFAULT_CELL_RADIUS;
		unordered_map<std::pair<int, int>, Cell*> ybuffer;
		for (auto& c : w.cells) {
			const auto& cpos = c->getPosition();
			// cellBL & cellUR are in sensors space
			std::pair<int, int> cellBottomLeft = {
			    floor((c->getPosition().x() - c->getBoundingBoxRadius()) / resolution),
			    floor((c->getPosition().z() - c->getBoundingBoxRadius()) / resolution)};
			std::pair<int, int> cellUpRight = {
			    ceil((c->getPosition().x() + c->getBoundingBoxRadius()) / resolution),
			    ceil((c->getPosition().z() + c->getBoundingBoxRadius()) / resolution)};
			for (int x = cellBottomLeft.first; x < cellUpRight.first; ++x) {
				for (int z = cellBottomLeft.second; z < cellUpRight.second; ++z) {
					std::pair<int, int> xz(x, z);
					if (cpos.y() > Config::EPSILON_GROUND &&
					    (!ybuffer.count(xz) || ybuffer[xz]->getPosition().y() < cpos.y()))
						ybuffer[xz] = c;
				}
			}
			c->setSensedNutrients(LIGHT, 0.0);
		}
		for (auto& c : ybuffer) {
			double lightIntensity = max(
			    0.0, Config::SUN_INTENSITY *
			             (min(1.0, c.second->getPosition().y() / Config::MAX_LIGHT_THRESHOLD)));
			c.second->setSensedNutrients(LIGHT, lightIntensity);
		}
	}

	void diffuseNutrients() {
		const auto dt = w.getDt();
		for (auto& c : w.cells) {
			double freeArea = c->getMembrane().getCurrentArea();
			for (const auto& con :
			     c->getMembrane().getCellCellConnectionManager().cellConnections) {
				freeArea -= con->area;
			}
			freeArea = std::max(0.0, freeArea);
			// light to cell
			if (c->sensedNutrients[LIGHT] > c->nutrientLevel[LIGHT]) {  // absorption only
				double deltaP =
				    max(0.0, c->nutrientLevel[LIGHT]) - max(0.0, c->sensedNutrients[LIGHT]);
				auto Qn = diffusedQtty(deltaP, freeArea, c->getMembrane().getBaseRadius(), dt);
				if (abs(Qn) > 10 || c->nutrientLevel[LIGHT] > 10) {
					std::cerr << "diffusing light to cell " << c->id
					          << ". c->lightLevel = " << c->nutrientLevel[LIGHT]
					          << ", freeArea = " << freeArea
					          << ", sensed = " << c->sensedNutrients[LIGHT]
					          << ", deltaP = " << deltaP << ", Qn = " << Qn << std::endl;
				}
				c->deltaNutrient(LIGHT, Qn);
			}
			// water to cell
			for (auto& wn : nutrientSources) {
				if (wn.content > 0) {
					auto pwater = computeNutrientIntensity(c->getPosition(), wn);
					if (pwater > c->nutrientLevel[WATER]) {  // absorption only
						double deltaP = max(0.0, c->nutrientLevel[WATER]) - max(0.0, pwater);
						auto Qn =
						    diffusedQtty(deltaP, freeArea, c->getMembrane().getBaseRadius(), dt);
						c->deltaNutrient(WATER, Qn);
						if (abs(Qn) > 10 || c->nutrientLevel[WATER] > 10) {
							std::cerr << "diffusing nutrients from grnd to cell " << c->id
							          << ". c->nutrientLevel = " << c->nutrientLevel[WATER]
							          << ", sensed = " << c->sensedNutrients[WATER]
							          << ", deltaP = " << deltaP << ", Qn = " << Qn << std::endl;
						}
						wn.content -= Qn;
					}
				}
			}
		}
		// cell to cell
		for (size_t n = 0; n < Config::NB_NUTRIENTS; ++n) {
			for (auto& con : w.cellCellConnections) {
				if (con.second->adhCoef > 0.1) {  // true connection
					double deltaP = max(0.0, con.first.second->nutrientLevel[n]) -
					                max(0.0, con.first.first->nutrientLevel[n]);
					auto Qn =
					    diffusedQtty(deltaP, con.second->adhArea, con.second->centersDist, dt);
					// if (abs(Qn) > 0) {
					// std::cerr << "diffusing nutrients " << n << " from cell"
					//<< con.second->cells.first->id << " to cell "
					//<< con.second->cells.second->id << ", deltaP = " << deltaP
					//<< ", Qn = " << Qn << std::endl;
					// std::cerr << " first cell level = " << con.first.first->nutrientLevel[n]
					//<< ", ";
					// std::cerr << " secnd cell level = " << con.first.second->nutrientLevel[n]
					//<< std::endl;
					con.first.second->deltaNutrient(n, Qn);
					con.first.first->deltaNutrient(n, -Qn);
					// std::cerr << "AFTERDIF:" << std::endl;
					// std::cerr << " first cell level = " << con.first.first->nutrientLevel[n]
					//<< ", ";
					// std::cerr << " secnd cell level = " << con.first.second->nutrientLevel[n]
					//<< std::endl;
					//}
				}
			}
		}
	}

	double computeNutrientIntensity(const MecaCell::Vec& p, const NutrientSource& n,
	                                bool sampling = false) const {
		if (p.y() > -Config::EPSILON_GROUND) return 0.0;
		double sqd = (p - n.pos).sqlength();
		// return max(0.0, 1.0 - (sqd / (n.sqradius * n.content / n.initialcontent)));
		// return max(0.0, 1.0 - (sqd / (n.sqradius * n.content / n.initialcontent)));
		double r = n.content / n.initialcontent;
		double coef = sampling ? Config::NUTRIENT_SAMPLING_COEF : 1.0;
		return max(0.0, n.content * (1.0 - (sqd / n.sqradius * r * coef)) * r);
	}
	double computeNutrientIntensity(const MecaCell::Vec& p) {
		double res = 0.0;
		for (const auto& n : nutrientSources) res += computeNutrientIntensity(p, n);
		return res;
	}

	double estimateFreeAreaSectionRatio(const MecaCell::Vec& dir, Cell* c) {
		double sumArea = 0.0;
		// for (const auto& con :
		// c->getMembrane().getCellCellConnectionManager().cellConnections) {
		// sumArea += con->area *
		// max(0.0, (c == con->cells.first ? con->normal : -con->normal).dot(dir));
		//}
		for (const auto& con :
		     c->getMembrane().getCellCellConnectionManager().cellConnections) {
			sumArea += con->area;
		}
		double baseArea = 0.8 * M_PI * pow(c->getBoundingBoxRadius(), 2);
		if (sumArea >= baseArea)
			return 0.0;
		else
			return (baseArea - sumArea) / baseArea;
	}

	void applyDrag(double visco, Cell* c) {
		// we apply drag using stokes' law but proportionnaly
		// to the free area section ratio in the direction of movement.
		c->receiveExternalForce(-6.0 * M_PI * visco * c->getBoundingBoxRadius() *
		                        c->getVelocity());
		// estimateFreeAreaSectionRatio(c->getVelocity().normalized(), c));
	}

	void applyGravity(double g) {
		for (auto& c : w.cells)
			if (c->getPosition().y() > 0)
				c->receiveExternalForce(MecaCell::Vec(0, -g * c->getMass(), 0));
	}

	void applyGroundReaction(Cell* c) {
		double area = std::pow(c->getBoundingBoxRadius(), 2) * M_PI;
		MecaCell::Vec f = c->getExternalForces() + c->getForce();
		double r = 1.0;  // estimateFreeAreaSectionRatio(f.normalized(), c);
		auto actualF = r * f;
		auto l = actualF.length();
		if (l <= r * Config::GROUND_REACTION * area) {
			c->setForce(MecaCell::Vec::zero());
		} else {
			c->receiveForce(-(actualF / l) * r * Config::GROUND_REACTION * area);
		}
	}

	void worldupdate() {
		w.prepareCellForNextUpdate();
		applyGravity(Config::GRAVITY);
		for (auto& c : w.cells) {
			if (c->getPosition().y() < 0) {
				applyGroundReaction(c);
			} else {
				applyDrag(Config::AIR_VISCOSITY, c);
			}
		}
		w.updateExistingCollisionsAndConnections();
		for (auto& c : w.cells) {
			if (c->getPosition().y() < 0)
				c->template updatePositionsAndOrientations<PosIntegrator>(Config::SIM_DT);
			else
				c->template updatePositionsAndOrientations<MecaCell::Euler>(Config::SIM_DT);
		}
		w.lookForNewCollisionsAndConnections();
		w.updateBehaviors();
		w.destroyDeadCells();
		w.frame++;
	}

	void loop() {
		simTime += Config::SIM_DT;
		updateCellsSensedNutrients();
		shineOn();
		diffuseNutrients();
		std::vector<std::array<std::pair<MecaCell::Vec, double>, Config::NB_MORPHOGENS>>
		    morphogens;
		cellgrid.clear();
		for (auto& c : w.cells) cellgrid.insertOnlyCenter(c);
		auto gridcontent = cellgrid.getContent();
		for (auto& gridcell : gridcontent) {
			std::array<std::pair<MecaCell::Vec, double>, Config::NB_MORPHOGENS> morphoCenters{};
			for (auto& c : gridcell.second) {
				for (auto i = 0u; i < Config::NB_MORPHOGENS; ++i) {
					morphoCenters[i].first += c->getPosition();
					morphoCenters[i].second += c->morphogensProduction[i];
				}
			}
			if (gridcell.second.size() > 0) {
				for (auto i = 0u; i < Config::NB_MORPHOGENS; ++i) {
					morphoCenters[i].first /= static_cast<double>(gridcell.second.size());
					morphoCenters[i].second /= static_cast<double>(gridcell.second.size());
				}
			}
			morphogens.push_back(morphoCenters);
		}
		for (auto& c : w.cells) c->updateInputs(morphogens, this);
		worldupdate();
	}

	void printState() {
		std::cerr << " ------------------------------ " << std::endl;
		std::cerr << " simTime = " << simTime << ", " << w.cells.size() << " cells"
		          << std::endl;
		for (auto& c : w.cells) {
			std::cerr << c->toString() << std::endl;
		}
	}

	World& getWorld() { return w; }

	bool finished() {
		if (w.cells.size() == 0) return true;
		if (w.cells.size() > maxCells) return true;
		if (simTime > simDuration) return true;
		return false;
	}
};
#endif
