#include "core/evaluators.hpp"
#include "core/typesconfig.hpp"
#include "external/cxxopts.hpp"
#include "external/gaga/gaga.hpp"
#include <sys/types.h>
#include <unistd.h>
using G = GRN<RealCoords>;
using MG = MGRN<MGClassic>;
using C = GRNPlantController<G>;
using MC = GRNPlantController<MG>;
template <typename Ctrl> using PC = PlantCell<Ctrl, MecaCell::VolumeMembrane>;

template <typename S, typename GA>
int launchGA(GA&& evo, const std::string& evaluatorName, const std::string& jobId, int argc, char** argv) {
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
	evo.setSaveFolder(std::string("evo_") + jobId + std::string("/"));
	evo.setVerbosity(1);
	evo.setPopSize(200);
	evo.setMutationProba(0.8);
	evo.setPopSaveInterval(100);
	evo.setGenSaveInterval(10);
	evo.setSaveIndStats(true);
	evo.setSaveGenStats(true);
	evo.initPopulation(
	    []() { return std::remove_reference<GA>::type::DNA_t::random(0, nullptr); });
	evo.setCrossoverProba(0.3);
	evo.setNbSavedElites(3);
	evo.step(700);
	return 0;
}

int main(int argc, char** argv) {
	std::string evaluatorName;
	std::string grnType;
	std::string jobId = std::to_string(getpid());
	try {
		cxxopts::Options options(argv[0]);
		options.add_options()("e,evaluator", "evaluator name",
		                      cxxopts::value<std::string>(evaluatorName));
		options.add_options()("g,grn", "grn type", cxxopts::value<std::string>(grnType));
		options.add_options()("j,job", "job id", cxxopts::value<std::string>(jobId));
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
		launchGA<Scenario<PC<C>>>(evo, evaluatorName, jobId, argc, argv);
	} else if (grnType == "mgrn") {
		GAGA::GA<GRNPlantController<MGRN<MGClassic>>> evo(argc, argv);
		launchGA<Scenario<PC<MC>>>(evo, evaluatorName, jobId, argc, argv);
	} else {
		std::cerr << "No valid grn type found, aborting." << std::endl;
		return 1;
	}
	return 0;
}
