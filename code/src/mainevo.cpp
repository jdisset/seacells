#include "core/evaluators.hpp"
#include "core/typesconfig.hpp"
#include "external/cxxopts.hpp"
#include "external/gaga/gaga.hpp"
using G = GRN<RealCoords>;
using MG = MGRN<MGClassic>;
using C = GRNPlantController<G>;
using MC = GRNPlantController<MG>;
template <typename Ctrl> using PC = PlantCell<Ctrl, MecaCell::VolumeMembrane>;

template <typename S, typename GA>
int launchGA(GA&& evo, const std::string& evaluatorName, int argc, char** argv) {
	if (evaluatorName == "survival_only") {
		evo.setEvaluator([=](auto& i) {
			SurvivalEvaluator<S> e(argc, argv);
			e(i);
		});
	} else if (evaluatorName == "survival_and_novelty") {
		evo.setEvaluator([=](auto& i) {
			SurvivalAndNoveltyEvaluator<S> e(argc, argv);
			e(i);
		});
		evo.enableNovelty();
		evo.setMinNoveltyForArchive(0.1);
	} else {
		std::cerr << "No valid evaluator found, aborting." << std::endl;
		return 1;
	}
	evo.setVerbosity(3);
	evo.setPopSize(200);
	evo.setMutationProba(0.9);
	evo.initPopulation(
	    []() { return std::remove_reference<GA>::type::DNA_t::random(0, nullptr); });
	// evo.setCrossoverProba(0.3);
	evo.setCrossoverProba(0.0);
	evo.setNbSavedElites(3);
	evo.step(400);
	return 0;
}

int main(int argc, char** argv) {
	std::string evaluatorName;
	std::string grnType;
	try {
		cxxopts::Options options(argv[0]);
		options.add_options()("e,evaluator", "evaluator name",
		                      cxxopts::value<std::string>(evaluatorName));
		options.add_options()("g,grn", "grn type", cxxopts::value<std::string>(grnType));
		options.parse(argc, argv);
	} catch (const cxxopts::OptionException& e) {
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	} catch (const std::bad_cast& e) {
		std::cout << "bad cast: " << e.what() << std::endl;
		exit(1);
	}

	if (grnType == "grn") {
		GAGA::GA<C> evo(argc, argv);
		launchGA<Scenario<PC<C>>>(evo, evaluatorName, argc, argv);
	} else if (grnType == "mgrn") {
		GAGA::GA<GRNPlantController<MGRN<MGClassic>>> evo(argc, argv);
		launchGA<Scenario<PC<MC>>>(evo, evaluatorName, argc, argv);
	} else {
		std::cerr << "No valid grn type found, aborting." << std::endl;
		return 1;
	}
	return 0;
}
