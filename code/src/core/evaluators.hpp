#include <string>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <numeric>
#include <deque>
#include <mecacell/mecacell.h>
#include "config.hpp"
#include "../external/grgen/common.h"
#include "capture.hpp"

template <class Scenario> struct ComplexMorphologyEvaluator {
	const std::string name = "complexMorpho";
	int argc;
	char **argv;

	ComplexMorphologyEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename Cell> double computeSphericity(const unordered_set<Cell *> clust) {
		MecaCell::Grid<Cell *> grid(MecaCell::DEFAULT_CELL_RADIUS / 3.0);
		for (auto &c : clust) grid.insert(c);
		auto s = grid.computeSphericity();
		return s;
	}

	std::vector<double> getFootprint(Scenario &sc) {
		if (sc.getWorld().cells.size() == 0) {
			return {{0, 0}};
		}
		auto clusters =
		    ClusterTools::getClusters(sc.getWorld().cells, sc.getWorld().cellCellConnections);
		size_t biggestClusterId = 0;
		size_t biggestClusterSize = clusters[biggestClusterId].size();
		for (size_t i = 1; i < clusters.size(); ++i) {
			if (clusters[i].size() > biggestClusterSize) {
				biggestClusterSize = clusters[i].size();
				biggestClusterId = i;
			}
		}
		return {{static_cast<double>(biggestClusterSize), 0.0}};

		if (sc.getWorld().cells.size() == 0) {
			return {{-1, -1, -1, -1}};
		}

		auto sphericity = computeSphericity(clusters[biggestClusterId]);
		vector<typename Scenario::CellType *> biggestClusterVec;
		biggestClusterVec.reserve(clusters[biggestClusterId].size());
		MecaCell::Vec centroid(0, 0, 0);
		for (auto &c : clusters[biggestClusterId]) {
			biggestClusterVec.push_back(c);
			centroid += c->getPosition();
		}
		if (biggestClusterSize > 0) {
			centroid /= static_cast<double>(biggestClusterSize);
		}

		auto longestAxis = getLongestDistance(biggestClusterVec);
		double longestProj2 = 0;
		double longestProj3 = 0;
		double longestDist = longestAxis.length();
		if (longestDist > 0) {
			MecaCell::Vec axis2 = longestAxis.ortho().normalized();
			MecaCell::Vec axis3 = longestAxis.cross(axis2).normalized();
			for (size_t i = 0; i < biggestClusterSize; ++i) {
				for (size_t j = i + 1; j < biggestClusterSize; ++j) {
					auto AB =
					    biggestClusterVec[i]->getPosition() - biggestClusterVec[j]->getPosition();
					auto proj2 = fabs(AB.dot(axis2));
					auto proj3 = fabs(AB.dot(axis3));
					if (longestProj2 < proj2) longestProj2 = proj2;
					if (longestProj3 < proj3) longestProj3 = proj3;
				}
			}
		}
		double normalizedNbCells = ClusterTools::sqAsympt(
		    biggestClusterSize, 0.002);  // counts less after 1000 cells.
		double ratioProj2 = 1.0;
		double ratioProj3 = 1.0;
		if (longestDist > 0) {
			ratioProj2 = longestProj2 / longestDist;
			ratioProj3 = longestProj3 / longestDist;
		}
		return {{normalizedNbCells, sphericity, ratioProj2, ratioProj3}};
	}

	template <typename Cell> MecaCell::Vec getLongestDistance(std::vector<Cell *> cells) {
		MecaCell::Vec longest(0, 0, 0);
		double sql = 0;
		for (size_t i = 0; i < cells.size(); ++i) {
			for (size_t j = 0; j < cells.size(); ++j) {
				auto axis = cells[i]->getPosition() - cells[j]->getPosition();
				auto l = axis.sqlength();
				if (l > sql) {
					sql = l;
					longest = axis;
				}
			}
		}
		return longest;
	}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		while (!sc.finished()) {
			sc.loop();
		}
		std::vector<std::vector<double>> footprints;
		sc.terminate();
		footprints.push_back(getFootprint(sc));
		ind.footprint = footprints;
	}
};

template <class Scenario> struct CaptureEvaluator {
	const std::string name = "Capture";
	int argc;
	char **argv;

	CaptureEvaluator(int c, char **v) : argc(c), argv(v) {}
	template <typename T>
	static std::vector<std::vector<double>> captMatrixTofootprint(const T &capture) {
		std::vector<std::vector<double>> res;
		res.reserve(capture.size());
		for (const auto &i : capture) {
			std::vector<double> line;
			line.reserve(i.size());
			for (auto j : i) {
				line.push_back(j);
			}
			res.push_back(line);
		}
		return res;
	}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		while (!sc.finished()) {
			sc.loop();
		}
		sc.terminate();
		const int W = 30;
		const int H = 30;
		std::ostringstream info;
		auto biggestCluster = ClusterTools::getBiggestCluster(
		    sc.getWorld().cells, sc.getWorld().cellCellConnections);
		auto capture = ClusterTools::getMatrixCapture<W, H>(
		    biggestCluster, true, 0);  // binarized img of the biggest cluster, auto scaled.
		typename Scenario::CellType stemCopy(ind.dna);
		auto nbC = biggestCluster.size();
		info << "  Total cells: " << sc.getWorld().cells.size() << ", biggest cluster:" << nbC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl
		     << ClusterTools::captMatrixToString(capture);
		ind.infos = info.str();
		ind.footprint = captMatrixTofootprint(capture);
	}
};

template <class Scenario> struct EnergyNoveltyEvaluator {
	const std::string name = "EnergyAndNovelty";
	int argc;
	char **argv;

	EnergyNoveltyEvaluator(int c, char **v) : argc(c), argv(v) {}
	template <typename T>
	static std::vector<std::vector<double>> captMatrixTofootprint(const T &capture) {
		std::vector<std::vector<double>> res;
		res.reserve(capture.size());
		for (const auto &i : capture) {
			std::vector<double> line;
			line.reserve(i.size());
			for (auto j : i) {
				line.push_back(j);
			}
			res.push_back(line);
		}
		return res;
	}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		while (!sc.finished()) {
			sc.loop();
		}
		sc.terminate();
		const int W = 30;
		const int H = 30;
		std::ostringstream info;
		auto biggestCluster = ClusterTools::getBiggestCluster(
		    sc.getWorld().cells, sc.getWorld().cellCellConnections);
		auto capture = ClusterTools::getMatrixCapture<W, H>(
		    biggestCluster, true, 0);  // binarized img of the biggest cluster, auto scaled.
		typename Scenario::CellType stemCopy(ind.dna);
		auto nbC = biggestCluster.size();
		info << "  Net energy: " << sc.plantEnergy
		     << ", total cells: " << sc.getWorld().cells.size() << ", biggest cluster:" << nbC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl
		     << ClusterTools::captMatrixToString(capture);
		ind.infos = info.str();
		ind.footprint = captMatrixTofootprint(capture);
		ind.fitnesses["Energy"] = sc.plantEnergy;
	}
};
template <class Scenario> struct EnergyEvaluator {
	const std::string name = "Energy";
	int argc;
	char **argv;

	EnergyEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		while (!sc.finished()) {
			sc.loop();
		}
		sc.terminate();
		const int W = 17;
		const int H = 17;
		std::ostringstream info;
		auto capture = ClusterTools::getMatrixCapture<W, H>(sc.getWorld().cells, false, 0);
		typename Scenario::CellType stemCopy(ind.dna);
		auto biggestCluster = ClusterTools::getBiggestCluster(
		    sc.getWorld().cells, sc.getWorld().cellCellConnections);
		auto nbC = biggestCluster.size();
		info << "  Net energy: " << sc.plantEnergy
		     << ", total cells: " << sc.getWorld().cells.size() << ", biggest cluster:" << nbC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl
		     << ClusterTools::captMatrixToString(capture);
		ind.infos = info.str();
		ind.fitnesses["Energy"] = sc.plantEnergy;
	}
};
template <class Scenario> struct SurvivalEvaluator {
	const std::string name = "Survival";
	int argc;
	char **argv;

	SurvivalEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		size_t maxC = 1;
		while (!sc.finished()) {
			sc.loop();
			if (sc.getWorld().cells.size() > maxC) maxC = sc.getWorld().cells.size();
		}
		sc.terminate();
		std::ostringstream info;
		typename Scenario::CellType stemCopy(ind.dna);
		info << "  Net energy: " << sc.plantEnergy << ", max cells: " << maxC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl;
		ind.infos = info.str();
		ind.fitnesses["Survival"] = sc.simTime;
	}
};

template <class Scenario> struct SurvivalNoveltyOnlyEvaluator {
	const std::string name = "SurvivalNoveltyOnly";
	int argc;
	char **argv;

	SurvivalNoveltyOnlyEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		size_t maxC = 1;
		while (!sc.finished()) {
			sc.loop();
			if (sc.getWorld().cells.size() > maxC) maxC = sc.getWorld().cells.size();
		}
		sc.terminate();
		std::ostringstream info;
		typename Scenario::CellType stemCopy(ind.dna);
		info << "  Survival time : " << sc.simTime << ", max cells: " << maxC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl;
		ind.infos = info.str();
		vector<vector<double>> fp;
		fp.push_back(vector<double>());
		fp[0].push_back(sc.simTime);
		ind.footprint = fp;
	}
};

template <class Scenario> struct SurvivalAndNoveltyEvaluator {
	const std::string name = "SurvivalAndNovelty";
	int argc;
	char **argv;

	SurvivalAndNoveltyEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename Individu> void operator()(Individu &ind) {
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		size_t maxC = 1;
		double minY = 0.0;
		while (!sc.finished()) {
			sc.loop();
			if (sc.getWorld().cells.size() > maxC) maxC = sc.getWorld().cells.size();
			if (sc.getWorld().cells.size() > 0) {
				for (auto &c : sc.getWorld().cells) {
					auto y = c->getPosition().y();
					if (y < minY && y > -1500) minY = y;
				}
			}
		}
		sc.terminate();
		std::ostringstream info;
		typename Scenario::CellType stemCopy(ind.dna);
		info << "  Survival time : " << sc.simTime << ", max cells: " << maxC
		     << ", max depth: " << minY << ", grn size: ["
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl;
		ind.infos = info.str();
		vector<vector<double>> fp;
		fp.push_back(vector<double>());
		fp[0].push_back(sc.simTime);
		fp[0].push_back(maxC);
		fp[0].push_back(minY / 40.0);
		ind.footprint = fp;
		ind.fitnesses["Survival"] = sc.simTime;
	}
};

template <class Scenario> struct SurvivalAndCaptureEvaluator {
	const std::string name = "SurvivalAndCapture";
	int argc;
	char **argv;

	SurvivalAndCaptureEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename T>
	static std::vector<double> captMatrixTofootprint(const T &capture) {
		// we flatten the image
		std::vector<double> res;
		res.reserve(capture.size() * capture.size());
		for (const auto &i : capture)
			for (auto j : i) res.push_back(j);
		return res;
	}

	template <typename Individu> void operator()(Individu &ind) {
		const std::array<double, 5> capturesTime = {{10.0, 20.0, 40.0, 60.0, 100.0}};
		const int W = 15;
		const int H = 15;
		const double maxH = 500.0;
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		ind.footprint = std::vector<vector<double>>();
		size_t maxC = 1;
		std::vector<std::string> capturesStr;
		while (!sc.finished()) {
			sc.loop();
			if (sc.getWorld().cells.size() > maxC) maxC = sc.getWorld().cells.size();
			while (ind.footprint.size() < capturesTime.size() &&
			       sc.simTime > capturesTime[ind.footprint.size()]) {
				auto capture =
				    ClusterTools::getMatrixCapture<W, H>(sc.getWorld().cells, false, maxH);
				ind.footprint.push_back(captMatrixTofootprint(capture));
				capturesStr.push_back(ClusterTools::captMatrixToString(capture));
			}
		}
		while (ind.footprint.size() < capturesTime.size()) {
			auto capture =
			    ClusterTools::getMatrixCapture<W, H>(sc.getWorld().cells, false, maxH);
			ind.footprint.push_back(captMatrixTofootprint(capture));
			capturesStr.push_back(ClusterTools::captMatrixToString(capture));
		}
		sc.terminate();
		std::ostringstream info;
		typename Scenario::CellType stemCopy(ind.dna);
		info << "  Survival time : " << sc.simTime << ", max cells: " << maxC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl;
		for (auto &cap : capturesStr) {
			info << cap << std::endl;
		}
		ind.infos = info.str();
		ind.fitnesses["Survival"] = sc.simTime;
	}
};

template <class Scenario> struct SurvivalAndMultiNoveltyEvaluator {
	const std::string name = "SurvivalAndMultiNovelty";
	int argc;
	char **argv;

	SurvivalAndMultiNoveltyEvaluator(int c, char **v) : argc(c), argv(v) {}

	template <typename Individu> void operator()(Individu &ind) {
		const std::array<double, 5> capturesTime = {{10.0, 20.0, 40.0, 60.0, 100.0}};
		Scenario sc;
		sc.setStemCell(new typename Scenario::CellType(ind.dna));
		sc.init(argc, argv);
		ind.footprint = std::vector<vector<double>>();
		size_t maxC = 1;
		while (!sc.finished()) {
			sc.loop();
			if (sc.getWorld().cells.size() > maxC) maxC = sc.getWorld().cells.size();
			while (ind.footprint.size() < capturesTime.size() &&
			       sc.simTime > capturesTime[ind.footprint.size()]) {
				vector<double> capture;
				capture.push_back(sc.getWorld().cells.size());

				if (sc.getWorld().cells.size() > 0) {
					double minY = 10000.0;
					for (auto &c : sc.getWorld().cells) {
						auto y = c->getPosition().y();
						if (y < minY && y > -1500) minY = y;
					}
					capture.push_back(minY / 40.0);
				} else
					capture.push_back(0.0);
				ind.footprint.push_back(capture);
			}
		}
		while (ind.footprint.size() < capturesTime.size()) {
			vector<double> capture;
			capture.push_back(-10);
			capture.push_back(0);
			ind.footprint.push_back(capture);
		}

		sc.terminate();
		std::ostringstream info;
		typename Scenario::CellType stemCopy(ind.dna);
		info << "  Survival time : " << sc.simTime << ", max cells: " << maxC
		     << ", grn size: [" << stemCopy.ctrl.grn.getProteinSize(ProteinType::input) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::regul) << ","
		     << stemCopy.ctrl.grn.getProteinSize(ProteinType::output) << "]" << endl;
		for (size_t i = 0; i < capturesTime.size(); ++i) {
			info << capturesTime[i] << " : " << ind.footprint[i][0] << ", "
			     << ind.footprint[i][1] << std::endl;
		}
		ind.infos = info.str();
		ind.fitnesses["Survival"] = sc.simTime;
	}
};
