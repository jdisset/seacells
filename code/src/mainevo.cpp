#include "external/gaga/gaga/gaga.hpp"
#include "external/cxxopts.hpp"
#include "core/evaluators.hpp"
#include "core/typesconfig.hpp"

template <typename GA> int launchGA(GA&& evo) {
	evo.setVerbosity(2);
	evo.setPopSize(200);
	evo.setNbGenerations(400);
	evo.setMutationProba(0.75);
	evo.setCrossoverProba(0.3);
	evo.setNbSavedElites(3);
	return evo.start();
}

int main(int argc, char** argv) {
	using scenario_t = TypesConfig::ScenarioType;
	using ctrl_t = TypesConfig::CtrlType;

	std::string evaluatorName;
	try {
		cxxopts::Options options(argv[0]);
		options.add_options()("e,evaluator", "evaluator name",
		                      cxxopts::value<std::string>(evaluatorName));
		options.parse(argc, argv);
	} catch (const cxxopts::OptionException& e) {
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	} catch (const std::bad_cast& e) {
		std::cout << "bad cast: " << e.what() << std::endl;
		exit(1);
	}
	if (evaluatorName == "survival")
		return launchGA(GAGA::GA<ctrl_t, SurvivalEvaluator<scenario_t>>(argc, argv));
	if (evaluatorName == "survival_novelty_only") {
		GAGA::GA<ctrl_t, SurvivalNoveltyOnlyEvaluator<scenario_t>> evo(argc, argv);
		evo.enableNovelty();
		evo.setMinNoveltyForArchive(0.1);
		return launchGA(evo);
	}
	if (evaluatorName == "survival_and_novelty") {
		GAGA::GA<ctrl_t, SurvivalAndNoveltyEvaluator<scenario_t>> evo(argc, argv);
		evo.enableNovelty();
		evo.setMinNoveltyForArchive(0.1);
		return launchGA(evo);
	}
	if (evaluatorName == "survival_and_capture") {
		GAGA::GA<ctrl_t, SurvivalAndCaptureEvaluator<scenario_t>> evo(argc, argv);
		evo.enableNovelty();
		evo.setMinNoveltyForArchive(3.0);
		return launchGA(evo);
	}
	if (evaluatorName == "survival_multinovelty") {
		GAGA::GA<ctrl_t, SurvivalAndMultiNoveltyEvaluator<scenario_t>> evo(argc, argv);
		evo.enableNovelty();
		evo.setMinNoveltyForArchive(1.0);
		return launchGA(evo);
	}

	std::cerr << "No valid evaluator found, aborting." << std::endl;
	return 0;
}
