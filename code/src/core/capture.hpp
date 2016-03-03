#ifndef CAPTURE_HPP
#define CAPTURE_HPP
#include <vector>
#include <string>
#include <deque>
#include <numeric>
#include <cmath>
#include <unordered_set>

struct ClusterTools {
	template <typename Cell, typename CC>
	static std::vector<std::unordered_set<Cell *>> getClusters(
	    const std::vector<Cell *> &cells, const CC &) {
		std::unordered_set<Cell *> clusterizedCells;
		std::vector<std::unordered_set<Cell *>> clusters;
		for (auto &c : cells) {
			if (!clusterizedCells.count(c)) {
				std::unordered_set<Cell *> clust;
				std::unordered_set<Cell *> toBeAdded = {c};
				while (toBeAdded.size()) {
					auto current = *toBeAdded.begin();
					clust.insert(current);
					toBeAdded.erase(toBeAdded.begin());
					for (auto &candidate : current->getConnectedCells()) {
						if (!clust.count(candidate)) {
							toBeAdded.insert(candidate);
						}
					}
				}
				for (auto &clustCell : clust) {
					clusterizedCells.insert(clustCell);
				}
				clusters.push_back(clust);
			}
		}
		return clusters;
	}

	static std::array<int, 3> HSVToRGB(float_t H, float_t S, float_t V) {
		// h =	[0, 360]; s, v = [0, 1]
		float_t C = V * S;  // Chroma
		float_t HPrime = fmod(H / 60.0, 6.0);
		float_t X = C * (1.0 - fabs(fmod(HPrime, 2.0) - 1.0));
		float_t M = V - C;
		float_t R, G, B;
		if (0 <= HPrime && HPrime < 1) {
			R = C;
			G = X;
			B = 0;
		} else if (1 <= HPrime && HPrime < 2) {
			R = X;
			G = C;
			B = 0;
		} else if (2 <= HPrime && HPrime < 3) {
			R = 0;
			G = C;
			B = X;
		} else if (3 <= HPrime && HPrime < 4) {
			R = 0;
			G = X;
			B = C;
		} else if (4 <= HPrime && HPrime < 5) {
			R = X;
			G = 0;
			B = C;
		} else if (5 <= HPrime && HPrime < 6) {
			R = C;
			G = 0;
			B = X;
		} else {
			R = 0;
			G = 0;
			B = 0;
		}
		R += M;
		G += M;
		B += M;
		return {{static_cast<int>(R * 255.0), static_cast<int>(G * 255.0),
		         static_cast<int>(B * 255.0)}};
	}
	static inline double sqAsympt(double x, double lambda) {
		auto a = lambda * lambda * x * x;
		return a / (a + 1.0);
	}
	template <typename Cell, typename CC>
	static vector<Cell *> getBiggestCluster(const vector<Cell *> &cells, const CC &cc) {
		vector<Cell *> biggestCluster;
		if (cells.size() > 0) {
			auto clusters = getClusters(cells, cc);
			size_t biggestClusterId = 0;
			size_t biggestClusterSize = clusters[biggestClusterId].size();
			for (size_t i = 1; i < clusters.size(); ++i) {
				if (clusters[i].size() > biggestClusterSize) {
					biggestClusterSize = clusters[i].size();
					biggestClusterId = i;
				}
			}
			biggestCluster.reserve(clusters[biggestClusterId].size());
			for (auto &c : clusters[biggestClusterId]) biggestCluster.push_back(c);
		}
		return biggestCluster;
	}

	template <typename Cell>
	static std::pair<MecaCell::Vec, MecaCell::Vec> fitPlane(vector<Cell *> cells) {
		// returns the best fitting plane {offset, normal} (min squared dist)
		size_t n = cells.size();
		if (n == 0) return {{0, 0, 0}, {0, 0, 1}};
		if (n == 1) return {cells[0]->getPosition(), {0, 0, 1}};
		if (n == 2) {
			auto AB = cells[1]->getPosition() - cells[0]->getPosition();
			auto l = AB.length();
			if (l > 0) AB /= l;
			return {cells[0]->getPosition() + AB * l * 0.5, AB.ortho()};
		}
		auto centroid = MecaCell::Vec::zero();
		for (auto &c : cells) centroid += c->getPosition();
		centroid /= static_cast<double>(n);

		double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;

		for (auto &c : cells) {
			auto p = c->getPosition() - centroid;
			xx += p.x() * p.x();
			xy += p.x() * p.y();
			xz += p.x() * p.z();
			yy += p.y() * p.y();
			yz += p.y() * p.z();
			zz += p.z() * p.z();
		}
		double detX = yy * zz - yz * yz;
		double detY = xx * zz - xz * xz;
		double detZ = xx * yy - xy * xy;

		double detMax = std::max(detX, std::max(detY, detZ));

		if (detMax <= 0.0) {  // not a plane
			auto AB = cells[1]->getPosition() - cells[0]->getPosition();
			auto l = AB.length();
			if (l > 0) AB /= l;
			return {cells[0]->getPosition() + AB * l * 0.5, AB.ortho()};
		}

		if (detMax == detX) {
			double a = (xz * yz - xy * zz) / detX;
			double b = (xy * yz - xz * yy) / detX;
			return {centroid, MecaCell::Vec(1.0, a, b).normalized()};
		} else if (detMax == detY) {
			double a = (yz * xz - xy * zz) / detY;
			double b = (xy * xz - yz * xx) / detY;
			return {centroid, MecaCell::Vec(a, 1.0, b).normalized()};
		} else {
			double a = (yz * xy - xz * yy) / detZ;
			double b = (xz * xy - yz * xx) / detZ;
			return {centroid, MecaCell::Vec(a, b, 1.0).normalized()};
		};
	}

	template <int FinalW = 30, int FinalH = 30, typename Cell>
	static std::array<std::array<double, FinalW>, FinalH> getMatrixCapture(
	    const vector<Cell *> &cells, bool binary = false, double maxWidth = 0.0) {
		// we oversample and downscale
		// we want W & H to be a multiple of FinalW & FinalH
		constexpr double overSamplingCoef = 3.0;
		constexpr int W = overSamplingCoef * FinalW;
		constexpr int H = overSamplingCoef * FinalH;

		std::array<std::array<double, H>, W> overSampledCapture = {};
		auto bestPlane = fitPlane(cells);
		const auto &centroid = bestPlane.first;
		const auto &zAxis = bestPlane.second;

		// finding main axis
		MecaCell::Vec xAxis(1, 0, 0);
		double sql = 0;
		for (size_t i = 0; i < cells.size(); ++i) {
			for (size_t j = 0; j < cells.size(); ++j) {
				auto axis = cells[j]->getPosition() - cells[i]->getPosition();
				auto projectedAxis = axis - zAxis * axis.dot(zAxis);
				auto l = projectedAxis.sqlength();
				if (l > sql) {
					sql = l;
					xAxis = projectedAxis;
				}
			}
		}
		double captWidth = maxWidth;
		if (maxWidth <= 0) {
			// auto scale, we asjust the so that all cells fit
			captWidth = sqrt(sql);
		}
		const double gridSize = (captWidth / static_cast<double>(W)) + 1.0;
		xAxis.normalize();
		auto yAxis = xAxis.cross(zAxis).normalized();
		// we align the longest axis along the bottomleft-upright diagonal
		// and we project each cell's bounding box in this plane
		for (auto &c : cells) {
			auto OC = c->getPosition() - centroid;
			double cellCenter_x = static_cast<double>(W / 2) * gridSize + OC.dot(xAxis);
			double cellCenter_y = static_cast<double>(H / 2) * gridSize + OC.dot(yAxis);
			int cellCenter_ix = floor(cellCenter_x / gridSize);
			int cellCenter_iy = floor(cellCenter_y / gridSize);
			double radius = c->getBoundingBoxRadius();
			int intRadius = ceil(radius / gridSize) + 1;
			for (int i = max(0, cellCenter_ix - intRadius);
			     i < min((int)W, cellCenter_ix + intRadius); ++i) {
				for (int j = max(0, cellCenter_iy - intRadius);
				     j < min((int)H, cellCenter_iy + intRadius); ++j) {
					double ipos = (static_cast<double>(i) + 0.5) * gridSize;
					double jpos = (static_cast<double>(j) + 0.5) * gridSize;
					double dist = sqrt(pow(ipos - cellCenter_x, 2) + pow(jpos - cellCenter_y, 2));
					double intensity = 1.0;
					if (dist > radius) intensity = 0.0;
					overSampledCapture[i][j] += intensity;
				}
			}
		}
		// downscaling using average
		std::array<std::array<double, FinalW>, FinalH> capture = {};
		for (int l = 0; l < FinalH; ++l) {
			for (int c = 0; c < FinalW; ++c) {
				double avg = 0;
				for (int x = c * overSamplingCoef; x < c * overSamplingCoef + overSamplingCoef;
				     ++x) {
					for (int y = l * overSamplingCoef; y < l * overSamplingCoef + overSamplingCoef;
					     ++y) {
						avg += overSampledCapture[x][y];
					}
				}
				avg /= (overSamplingCoef * overSamplingCoef);
				if (binary)
					capture[l][c] = avg > 0.5 ? 1 : 0;  // binarisation
				else
					capture[l][c] = avg;
			}
		}
		return capture;
	}
	template <typename T>
	static std::string captMatrixToString(const T &capture, bool trueColorEnabled = false) {
		const std::array<std::string, 8> pixels = {
		    {"⬣  ", "◉  ", "✹  ", "✷  ", "✺  ", "✳  ", "○  ", "⠂  "}};
		std::ostringstream output;
		output << "┏";

		for (size_t i = 0; i < capture[0].size(); ++i) {
			output << "---";
		}
		output << "┓" << std::endl;
		for (auto &i : capture) {
			std::ostringstream line;
			for (auto j : i) {
				if (j > 0.0) {
					int intensity =
					    7 - min(static_cast<int>(j * 5.0), 7);  // ouais c'est à l'envers...
					if (trueColorEnabled) {
						auto color = HSVToRGB(fabs(130.0 - j * 5.0), 0.75, 0.30 + min(j, 5.0) * 0.1);
						line << "\x1b[38;2;" << color[0] << ";" << color[1] << ";" << color[2] << "m"
						     << pixels[intensity] << "\x1b[0m";
					} else {
						line << pixels[intensity];
					}
				} else {
					line << "   ";
				}
			}
			output << "┇" << line.str() << "┇" << std::endl;
		}
		output << "┗";
		for (size_t i = 0; i < capture[0].size(); ++i) {
			output << "---";
		}
		output << "┛" << std::endl;
		return output.str();
	}
};

#endif
