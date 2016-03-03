#ifndef PLANTCELL_HPP
#define PLANTCELL_HPP
#include "config.hpp"
#include "../external/grgen/common.h"
#include <mecacell/mecacell.h>
#include <random>

enum class CycleStep { quiescent, growing };
template <typename Controller, template <class> class Membrane = MecaCell::VolumeMembrane>
class PlantCell
    : public MecaCell::ConnectableCell<PlantCell<Controller, Membrane>, Membrane> {
 public:
	using Base = MecaCell::ConnectableCell<PlantCell<Controller, Membrane>, Membrane>;
	using Vec = MecaCell::Vec;
	using CtrlType = Controller;
	using morphogrid =
	    std::vector<std::array<std::pair<MecaCell::Vec, double>, Config::NB_MORPHOGENS>>;

	std::normal_distribution<double> growthDistribution;
	std::mt19937 internalRand;
	std::array<double, Config::NB_MORPHOGENS> morphogensProduction{};
	std::array<double, Config::NB_MORPHOGENS> sensedMorphogens{};
	double age = 0.0;
	std::array<double, Config::NB_NUTRIENTS> sensedNutrients{};
	// sensed in the environment
	std::array<double, Config::NB_NUTRIENTS> nutrientLevel{};  // actual amount
	int needToComputeGradient = -1;
	double morphoUpdateDt = 0.0;
	CycleStep currentStep = CycleStep::quiescent;
	Controller ctrl;
	MecaCell::Vec divisionDirection{0, 0, 0};
	unordered_set<PlantCell*> trulyConnectedCells{};

	PlantCell(const Vec& p) : Base(p), ctrl(Controller::random(0, 0)) {
		init();
		nutrientLevel.fill(0.2);
	}

	PlantCell(const PlantCell& c, const Vec& p = Vec(0, 0, 0))
	    : Base(c, p),
	      morphogensProduction(c.morphogensProduction),
	      sensedMorphogens(c.sensedMorphogens),
	      age(0.0),
	      ctrl(c.ctrl) {
		init();
		for (size_t i = 0; i < nutrientLevel.size(); ++i) {
			nutrientLevel[i] = c.nutrientLevel[i];
			sensedNutrients[i] = c.sensedNutrients[i];
		}
	}

	PlantCell(const Controller& ct, const Vec& p = Vec(0, 0, 0)) : Base(p), ctrl(ct) {
		init();
		nutrientLevel.fill(0.2);
	}

	void init() {
		age = 0.0;
		sensedNutrients.fill(0.0);
		needToComputeGradient = -1;
		morphoUpdateDt = 0.0;
		currentStep = CycleStep::quiescent;
		divisionDirection = MecaCell::Vec(0, 0, 0);

		internalRand =
		    std::mt19937(static_cast<int>(this->getPosition().x() * this->getPosition().y() +
		                                  this->getPosition().z() + nutrientLevel[0] * 10.0));
		growthDistribution = std::normal_distribution<double>(
		    Config::CELL_GROWTH_SPEED * Config::SIM_DT,
		    Config::CELL_GROWTH_SPEED * Config::SIM_DT * 0.5);
	}

	void updateMorphogensProduction() {
		auto on = ctrl.getOutput("on");
		for (auto i = 0u; i < Config::NB_MORPHOGENS; ++i) {
			auto o = ctrl.getOutput(std::string("o") + std::to_string(i));
			morphogensProduction[i] = (on > 0.0 && o > 0.0) ? o / (o + on) : 0.0;
		}
	}

	double getAdhesionWith(PlantCell* c, MecaCell::Vec) const {
		if (Config::ENABLE_SOLIDIFY) {
			return (trulyConnectedCells.count(c) ||
			        ctrl.getOutput("s") < ctrl.getOutput("st")) ?
			           1.0 :
			           0.0;
		} else
			return 1.0;
	}

	double computeMorphogenIntensity(size_t i, const MecaCell::Vec& P,
	                                 const morphogrid& mg) {
		double sm = 0.0;
		for (const auto& c : mg) {
			auto sql = (c[i].first - P).sqlength() / Config::morphoDiffusionCoefs[i];
			sm += c[i].second / (sql + 1.0);
		}
		return sm;
	}

	void setSensedNutrients(size_t n, double cn) { sensedNutrients[n] = cn; }

	void deltaNutrient(size_t n, double amount) { nutrientLevel[n] += amount; }

	template <typename Sc> void updateInputs(const morphogrid& mg, const Sc* scenar) {
		if (morphoUpdateDt >= Config::MORPHOGEN_UPDATE_INTERVAL) morphoUpdateDt = 0.0;
		if (morphoUpdateDt == 0.0) {
			for (auto i = 0u; i < Config::NB_MORPHOGENS; ++i) {
				sensedMorphogens[i] = computeMorphogenIntensity(i, this->getPosition(), mg);
				ctrl.setInput(std::string("c") + std::to_string(i), sensedMorphogens[i]);
			}
		}
		for (auto i = 0u; i < Config::NB_NUTRIENTS; ++i) {
			ctrl.setInput(std::string("n") + std::to_string(i), nutrientLevel[i]);
			ctrl.setInput(std::string("cn") + std::to_string(i), sensedNutrients[i]);
		}
		auto normalizedAge = (0.05 * age) / (0.05 * age + 1.0);
		ctrl.setInput("t", normalizedAge);
		ctrl.setInput("bias", 1.0);
		ctrl.setInput("p", this->getNormalizedPressure());

		if (needToComputeGradient >= 0) {
			if (needToComputeGradient == Config::NB_MORPHOGENS) {
				divisionDirection = computeNutrientGradient(scenar);
			} else {
				auto gradient = computeMorphogenGradient(needToComputeGradient, mg);
				if (gradient.sqlength() > 0 &&
				    ctrl.getOutput("pd") >
				        ctrl.getOutput(std::string("d") +
				                       std::to_string(needToComputeGradient))) {
					// orthogonal division
					auto grad0 = computeMorphogenGradient(0, mg);
					if (grad0.sqlength() > 0 && needToComputeGradient > 0 &&
					    abs(gradient.dot(grad0)) < 0.99999999999) {
						divisionDirection = gradient.cross(grad0);
						if (divisionDirection.sqlength() > 0) divisionDirection.normalize();
					} else {
						std::normal_distribution<float_t> nDist(0.0, 1.0);
						MecaCell::Vec rdm(nDist(internalRand), nDist(internalRand),
						                  nDist(internalRand));
						divisionDirection = gradient.cross(rdm.normalized()).normalized();
					}
				} else {
					divisionDirection = gradient;
				}
			}
			needToComputeGradient = -1;
		}
	}

	void updateTrulyConnectedCells() {
		trulyConnectedCells = std::unordered_set<PlantCell*>();
		for (auto& con : this->membrane.getCellCellConnectionManager().cellConnections) {
			if (con->adhCoef > 0.1) {
				if (this == con->cells.first)
					trulyConnectedCells.insert(con->cells.second);
				else
					trulyConnectedCells.insert(con->cells.first);
			}
		}
	}

	PlantCell* updateBehavior(double dt) {
		for (auto i = 0u; i < nutrientLevel.size(); ++i) {
			if (nutrientLevel[i] < 0) {
				this->die();
				return nullptr;
			}
		}
		age += dt;
		morphoUpdateDt += dt;
		updateTrulyConnectedCells();
		ctrl.update();
		updateMorphogensProduction();
		this->setColorHSV(360.0 + 50.0 - nutrientLevel[WATER] * 200.0, 0.85,
		                  0.7 + 0.2 * sensedNutrients[LIGHT]);
		// this->color = {{sensedMorphogens[0] * 0.2, sensedMorphogens[1] * 0.2,
		// sensedMorphogens[2] * 0.2}};
		// this->color = {
		//{sensedNutrients[0], sensedNutrients[1], nutrientLevel[0] + nutrientLevel[1]}};
		if (currentStep == CycleStep::growing) {
			if (this->getRelativeVolume() >= 2.0) {
				if (divisionDirection.sqlength() == 0.0) {
					// no morpho gradient, random direction
					std::normal_distribution<float_t> nDist(0.0, 1.0);
					divisionDirection.coords = {
					    {nDist(internalRand), nDist(internalRand), nDist(internalRand)}};
					divisionDirection.normalize();
				}
				this->setMass(this->getBaseMass());
				this->membrane.division();
				this->currentStep = CycleStep::quiescent;
				for (auto& n : nutrientLevel) {
					n *= 0.5;  // dividing nutrients by 2
				}
				auto decalage =
				    divisionDirection * this->getMembraneDistance(divisionDirection) * 1.1;
				this->setPosition(this->getPosition() - (decalage * 0.1));
				auto nc = new PlantCell(*this, decalage);
				return nc;
			} else {
				this->grow(std::max(growthDistribution(internalRand),
				                    0.0));  // add a bit of randomness for asynchronicity
				for (auto& n : nutrientLevel) n -= Config::DIVISION_NRJ_CONSUMPTION * dt;
			}
		} else {
			currentStep = CycleStep::quiescent;
			size_t idStrongestDivGradient = 0;
			double maxDivOut = ctrl.getOutput("d0");
			for (auto i = 1u; i <= Config::NB_MORPHOGENS; ++i) {
				double concentration = 0.0;
				if (i == Config::NB_MORPHOGENS)
					concentration = ctrl.getOutput(std::string("dn"));
				else
					concentration = ctrl.getOutput(std::string("d") + std::to_string(i));
				if (concentration > maxDivOut) {
					maxDivOut = concentration;
					idStrongestDivGradient = i;
				}
			}
			double apop = ctrl.getOutput("a");
			double quiesc = ctrl.getOutput("q");
			if (maxDivOut > quiesc && maxDivOut > apop) {
				currentStep = CycleStep::growing;
				needToComputeGradient = idStrongestDivGradient;
			} else if (apop > quiesc) {
				this->die();
			}
			for (auto& n : nutrientLevel) n -= Config::NORMAL_NRJ_CONSUMPTION * dt;
		}
		return nullptr;
	}

	template <typename Sc> MecaCell::Vec computeNutrientGradient(const Sc* scenar) {
		using V = MecaCell::Vec;
		const auto d = Config::NUTRIENT_SAMPLING_DIST;
		const auto p = this->getPosition();
		MecaCell::Vec res(0, 0, 0);
		for (auto& n : scenar->nutrientSources) {
			res += MecaCell::Vec(scenar->computeNutrientIntensity(p + V(d, 0, 0), n, true) -
			                         scenar->computeNutrientIntensity(p - V(d, 0, 0), n, true),
			                     scenar->computeNutrientIntensity(p + V(0, d, 0), n, true) -
			                         scenar->computeNutrientIntensity(p - V(0, d, 0), n, true),
			                     scenar->computeNutrientIntensity(p + V(0, 0, d), n, true) -
			                         scenar->computeNutrientIntensity(p - V(0, 0, d), n, true));
		}
		res /= static_cast<double>(scenar->nutrientSources.size());
		if (res.sqlength() > 0)
			return -res.normalized();
		else
			return res;
	}

	MecaCell::Vec computeMorphogenGradient(size_t i, const morphogrid& mg) {
		using V = MecaCell::Vec;
		const auto d = Config::MORPHO_SAMPLING_DIST;
		const auto p = this->getPosition();
		MecaCell::Vec res(computeMorphogenIntensity(i, p + V(d, 0, 0), mg) -
		                      computeMorphogenIntensity(i, p - V(d, 0, 0), mg),
		                  computeMorphogenIntensity(i, p + V(0, d, 0), mg) -
		                      computeMorphogenIntensity(i, p - V(0, d, 0), mg),
		                  computeMorphogenIntensity(i, p + V(0, 0, d), mg) -
		                      computeMorphogenIntensity(i, p - V(0, 0, d), mg));
		if (res.sqlength() > 0)
			return res.normalized();
		else
			return res;
	}
};

#endif
