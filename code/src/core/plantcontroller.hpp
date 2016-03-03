#ifndef PLANTCONTROLLER_HPP
#define PLANTCONTROLLER_HPP
#include <string>
#include <sstream>
#include "config.hpp"
#include "../external/grgen/common.h"

template <typename GRN> struct GRNPlantController {
	static constexpr unsigned int nbMorphogens = Config::NB_MORPHOGENS;

	GRN grn;

	GRNPlantController() { reset(); }
	GRNPlantController(const GRN &g) : grn(g) { reset(); }
	GRNPlantController(const std::string &s) : grn(s) { reset(); }
	GRNPlantController(const GRNPlantController &other) : grn(other.grn) {}
	GRNPlantController &operator=(const GRNPlantController &other) {
		if (this != &other) {
			grn = other.grn;
			reset();
		}
		return *this;
	}

	void update() { grn.step(Config::GRN_STEPS_PER_UPDATE); }
	// GA specific methods
	GRNPlantController crossover(const GRNPlantController &other) {
		GRN g = grn.crossover(other.grn);
		GRNPlantController res(g);
		return res;
	}
	void mutate() { grn.mutate(); }
	void reset() { grn.reset(); }
	void setInput(const std::string &input, double val) {
		grn.setProteinConcentration(input, ProteinType::input, val);
	}
	double getOutput(const std::string &output) const {
		auto r = grn.getProteinConcentration(output, ProteinType::output);
		return r;
	}
	std::string toJSON() const { return grn.toJSON(); }

	static GRNPlantController random(int, char **) {
		// inputs
		GRN g;
		for (auto i = 0u; i < nbMorphogens; ++i)
			g.addRandomProtein(ProteinType::input, std::string("c") + std::to_string(i));
		for (auto i = 0u; i < Config::NB_NUTRIENTS; ++i) {
			g.addRandomProtein(ProteinType::input, std::string("n") + std::to_string(i));
			g.addRandomProtein(ProteinType::input, std::string("cn") + std::to_string(i));
		}
		g.addRandomProtein(ProteinType::input, std::string("t"));
		g.addRandomProtein(ProteinType::input, std::string("p"));
		g.addRandomProtein(ProteinType::input, std::string("bias"));
		// outputs
		for (auto i = 0u; i < nbMorphogens; ++i)
			g.addRandomProtein(ProteinType::output, std::string("o") + std::to_string(i));
		g.addRandomProtein(ProteinType::output, std::string("on"));
		for (auto i = 0u; i < nbMorphogens; ++i)
			g.addRandomProtein(ProteinType::output, std::string("d") + std::to_string(i));
		g.addRandomProtein(ProteinType::output, std::string("dn"));
		g.addRandomProtein(ProteinType::output, std::string("a"));
		g.addRandomProtein(ProteinType::output, std::string("q"));
		g.addRandomProtein(ProteinType::output, std::string("s"));
		g.addRandomProtein(ProteinType::output, std::string("st"));
		g.addRandomProtein(ProteinType::output, std::string("pd"));
		// reguls
		g.randomReguls(Config::INITIAL_NB_REGULS);
		g.randomParams();
		return GRNPlantController(g);
	}
};
#endif
